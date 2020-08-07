<script>
	import { fade } from 'svelte/transition'
  import { onMount } from 'svelte'
  import * as JsSearch from 'js-search'
  import mdIt from 'markdown-it'

  import { hash } from './stores'
  import Masonry from './Masonry'
  
  const base_url = window.location.origin

  // https://stackoverflow.com/questions/3426404/create-a-hexadecimal-colour-based-on-a-string-with-javascript
  function hashCode(str) {
      var hash = 0;
      for (var i = 0; i < str.length; i++) {
        hash = str.charCodeAt(i) + ((hash << 5) - hash);
      }
      return hash;
  } 

  function intToRGB(i){
      var c = (i & 0x00FFFFFF)
          .toString(16)
          .toUpperCase();
      return '#' + ("00000".substring(0, 6 - c.length) + c);
  }

  function localize_appyter_image({ name, image }) {
    if (/^https?:\/\//.exec(image) !== null) {
      return image
    } else {
      return `${base_url}/${name}/static/${image}`
    }
  }

  // store appyters as list and lookup table based on name slugs
  let appyterList = require('./appyters.json')
  // assemble appyter lookup table
  let appyterLookup = {}
  let globalTags = {}
  for (const appyter of appyterList) {
    let {name, description, long_description, authors, tags, ..._} = appyter
    const md = mdIt()
    const normalizeLink = md.normalizeLink
    md.normalizeLink = function (url) {
      if (/^https?:\/\//.exec(url) !== null) {
        return normalizeLink(url)
      } else if (/^\.\//.exec(url) !== null) {
        return normalizeLink(`${base_url}/${name}/${url.slice(2)}`)
      } else {
        return normalizeLink(`${base_url}/${name}/${url}`)
      }
    }
    const authors_flat = authors.map(({ name, email }) => `${name || ''} (${email || ''})`).join(', ')
    const description_html = md.render(description || '')
    // A bit roundabout but seemingly the easiest way to add img-fluid class to all markdown-rendered img tags
    let long_description_html = document.createElement('div')
    long_description_html.innerHTML = md.render((long_description || description).split('\n').slice(1).join('\n'))
    for (const img of long_description_html.querySelectorAll('img')) {
      img.classList.add('img-fluid')
    }
    long_description_html = long_description_html.innerHTML
    const color = intToRGB(hashCode(name))
    // modify appyters in-place
    Object.assign(appyter, {
      authors_flat,
      description_html,
      long_description_html,
      color,
    })
    // save tags to global list
    for (const tag of tags) {
      globalTags[tag] = (globalTags[tag] || 0) + 1
    }
    // save a reference in the lookup table
    appyterLookup[name] = appyter
  }
  for (const tag in globalTags) {
    if (globalTags[tag] <= 1) delete globalTags[tag]
  }
  appyterList = appyterList.filter(appyter => appyter.public !== false)

  // index documents for search
  const search = new JsSearch.Search('name')
  search.addIndex('name')
  search.addIndex('title')
  search.addIndex('authors_flat')
  search.addIndex('description')
  search.addIndex('long_description')
  search.addIndex('license')
  search.addIndex('tags')
  search.addIndex('url')
  search.addDocuments(appyterList)

  // get appyter hits
  async function get_pagehits() {
    const response = await fetch(
      `${base_url}/postgrest/pagehits?url=like.${encodeURIComponent(`${base_url}%`)}`
    )
    const pagehits = await response.json()
    for (const {url, hits} of pagehits) {
      if (url.startsWith(base_url + '/#')) {
        const appyter_name = url.slice(base_url.length + 2)
        if (appyterLookup[appyter_name] !== undefined) {
          Object.assign(appyterLookup[appyter_name], { views: hits })
        }
      } else if (url.startsWith(base_url + '/')) {
        if (url.endsWith('#view')) {
          const appyter_name = url.slice(base_url.length + 1, url.length - '#view'.length)
          if (appyterLookup[appyter_name] !== undefined) {
            Object.assign(appyterLookup[appyter_name], { persistent_views: hits })
          }
        } else if (url.endsWith('#execute')) {
          const appyter_name = url.slice(base_url.length + 1, url.length - '#execute'.length)
          if (appyterLookup[appyter_name] !== undefined) {
            Object.assign(appyterLookup[appyter_name], { runs: hits })
          }
        } else {
          const appyter_name = url.slice(base_url.length + 1)
          if (appyterLookup[appyter_name] !== undefined) {
            Object.assign(appyterLookup[appyter_name], { form_views: hits })
          }
        }
      }
    }
    appyterList.sort((a, b) => (b.views||0) - (a.views||0))
    appyterList = appyterList
  }

  async function pagehit(appyter) {
    await fetch(`${base_url}/postgrest/rpc/pagehit`, {
      method: 'post',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        'pageurl': appyter === undefined ? base_url : `${base_url}/#${appyter.name}`,
      }),
    })
  }

  function set_gte(A, B) {
    for (const b of B) {
      if (!A.has(b)) return false
    }
    return true
  }

  // facilitate search
  const searchAppyters = (appyterList, search, searchString, tags) => {
    let rs
    if (searchString !== '' && searchString !== undefined) {
      rs = search.search(searchString)
    } else {
      rs = appyterList
    }
    if (tags !== '' && tags !== undefined) {
      tags = new Set(tags.split(';'))
      rs = rs.filter(appyter => set_gte(new Set(appyter.tags), tags))
    }
    return rs
  }

  // sync appyter variable and url hash
  let appyter
  let lastPath
  $: {
    const curPath = $hash.path.slice(1)
    if (curPath !== lastPath) {
      appyter = appyterLookup[$hash.path.slice(1)]
      pagehit(appyter)
      lastPath = curPath
    }
    window.scrollTo(0,0)
  }

  // things to do on window load
  onMount(() => {
    get_pagehits()
  })
</script>

<style>
/* sticky footer */
:global(body) {
  display: flex;
  min-height: 100vh;
  flex-direction: column;
}

.flex-grow {
  flex: 1 0 auto;
}

.link-unstyled {
  color: inherit;
  text-decoration: inherit;
}

/* card shadow on hover */
.card {
  border-radius: 6px;
  box-shadow: 0 2px 4px 0 rgba(0,0,0,0.1);
  transition-timing-function: ease-in;
  transition: box-shadow 1s, transform 1s;
}
.card:hover {
  box-shadow: 0 8px 16px 0 rgba(0,0,0,0.2);
  transition-timing-function: ease-out;
  transition: box-shadow 1s;
}

:global(body) {
  background-color: #f5f5f5;
}
</style>

<div class="row m-0 p-0">
  <div class="col-sm-12 offset-md-3 col-md-6 offset-lg-4 col-lg-4 text-center">
    <h1>
      <a href="#/">
        <img
          src="{require('./images/appyters_logo.svg')}"
          class="img-fluid w-100 p-2"
          alt="Appyters"
        />
      </a>
    </h1>
    <h3 class="card-subtitle mb-2 text-muted">A catalog of appyter notebooks</h3>
  </div>
  <div class="col-sm-12 p-0">
    <hr />
  </div>
  {#if appyter === undefined}
    <div class="container">
      <div class="offset-sm-2 col-sm-8 text-center mb-3">
        <div class="row">
          <div class="col-sm-12">
            <input
              type="text"
              class="form-control"
              placeholder="Search appyters..."
              aria-label="Search appyters"
              bind:value={$hash.params.q}
            />
          </div>
          {#if $hash.params.tags}
            <div class="col-sm-12 mb-2">
              {#each $hash.params.tags.split(';') as tag}
                <a
                  href="javascript:"
                  class="badge badge-success m-1 p-1 text-white"
                  on:click={() => $hash.params.tags = $hash.params.tags.split(';').filter(t => t !== tag).join(';')}
                >{tag}</a>
              {/each}
            </div>
          {/if}
          <div class="col-sm-12">
            {#each Object.keys(globalTags).filter(tag => ($hash.params.tags || '').split(';').indexOf(tag) === -1) as tag}
              <a
                href="javascript:"
                class="badge badge-primary m-1 p-1 text-white"
                on:click={() => $hash.params.tags = [...($hash.params.tags || '').split(';').filter(t => t !== '' && t !== tag), tag].join(';')}
              >{tag}</a>
            {/each}
          </div>
        </div>
      </div>
    </div>
  {/if}
</div>
<div class="container content flex-grow">
  {#if appyter === undefined}
    <Masonry>
      {#each searchAppyters(appyterList, search, $hash.params.q, $hash.params.tags) as appyter}
        <div
          in:fade="{{duration: 250}}"
          out:fade="{{duration: 250}}"
        >
          <a
            href="#/{appyter.name}"
            class="link-unstyled"
          >
            <div class="card">
              <div
                class="card-img-top"
                style={[
                  `background-color: ${appyter.color}`,
                  appyter.image !== undefined ? (
                    `background-image: url('${localize_appyter_image(appyter)}')`
                  ) : undefined,
                  `background-repeat: no-repeat`,
                  `background-size: cover`,
                  `background-position: center`,
                  `width: 100%`,
                  `padding-top: 56.25%`, // 720 / 1280 = 0.5625, this preserves aspect ratio of div
                ].filter((el) => el !== undefined).join('; ')}
              ></div>
              <div class="card-body">
                <h3 class="card-title">{appyter.title}</h3>
                <div class="d-flex flex-row flex-nowrap pt-1 pb-2">
                  <div class="d-flex flex-column pr-2">
                  {#if appyter.views }
                    <div class="text-grey text-nowrap">
                      Views: {appyter.views}
                    </div>
                  {/if}
                  {#if appyter.runs }
                    <div class="text-grey text-nowrap">
                    Runs: {appyter.runs}
                    </div>
                  {/if}
                  </div>
                  <div class="d-flex flex-column pl-2">
                    {#if appyter.form_views }
                      <div class="text-grey text-nowrap">
                      Starts: {appyter.form_views}
                      </div>
                    {/if}
                    {#if appyter.persistent_views }
                      <div class="text-grey text-nowrap">
                      Retrievals: {appyter.persistent_views}
                      </div>
                    {/if}
                  </div>
                </div>
                <p class="card-text">{@html appyter.description_html}</p>
                <div class="pb-4 d-flex flex-row flex-wrap">
                  <span class="badge badge-success m-1 p-1">v{appyter.version}</span>
                  <a
                    class="badge badge-secondary m-1 p-1"
                    href="javascript:"
                    on:click={() => $hash.params.license = appyter.license}
                  >{appyter.license}</a>
                  {#each appyter.tags as tag}
                    <a
                      class="badge badge-primary m-1 p-1 text-white"
                      href="javascript:"
                      on:click={() => $hash.params.tags = [...($hash.params.tags || '').split(';').filter(t => t !== '' && t !== tag), tag].join(';')}
                    >{tag}</a>
                  {/each}
                </div>
                <button class="btn btn-primary btn-sm">
                  Select
                </button>
              </div>
            </div>
          </a>
        </div>
      {/each}
    </Masonry>
  {:else}
    <div class="row">
      <div class="col-sm-12 col-md-6 col-lg-4 align-self-center">
        <div 
          class="card-img-top"
          style={[
            `background-color: ${appyter.color}`,
            appyter.image !== undefined ? (
              `background-image: url('${localize_appyter_image(appyter)}')`
            ) : undefined,
            `background-repeat: no-repeat`,
            `background-size: cover`,
            `background-position: center`,
            `width: 100%`,
            `padding-top: 56.25%`, // 720 / 1280 = 0.5625, this preserves aspect ratio of div
          ].filter((el) => el !== undefined).join('; ')}
        >
        </div>
      </div>
      <div class="col-sm-12 col-md-6 col-lg-8">
        <h2>{appyter.title}</h2>
        <span class="card-subtitle mb-2 text-muted d-flex flex-row flex-wrap">
          <span class="badge badge-success m-1">v{appyter.version}</span>
          <span class="badge badge-secondary m-1">{appyter.license}</span>
          {#each appyter.tags as tag}
            <span class="badge badge-primary m-1">{tag}</span>
          {/each}
        </span>
        {#if appyter.url !== undefined}
          <p><a href="{appyter.url}">{appyter.url}</a></p>
        {/if}
        <p>
          <b>First Published:</b> {(new Date(appyter.creation_timestamp)).toDateString()}<br />
          <b>Last Updated:</b> {(new Date(appyter.update_timestamp)).toDateString()}<br />
        </p>
        <p>
          <b>Authors:</b><br />
          {#each appyter.authors as author}
            <span>{author.name} &lt;<a href="mailto:{author.email}">{author.email}</a>&gt;</span><br />
          {/each}
        </p>
      </div>
      <div class="col-sm-12">
        {@html appyter.long_description_html}
      </div>
      <div class="col-sm-12">
        <a href="{base_url}/{appyter.name}/" class="btn btn-primary">Start Appyter</a>
      </div>
    </div>
  {/if}
</div>
<div class="footer" style="margin-top: 25px;">
  <hr />
  <div class="row justify-content-center align-content-center" style="margin: 0">
    <div class="col-md-3 col-sm-12 text-center">
      <p>
        <a style="color: #555;" href="mailto:avi.maayan@mssm.edu">Contact Us</a><br />
        <a style="color: #555;" href="https://github.com/MaayanLab/appyter-catalog/blob/master/LICENSE">Usage License</a><br />
        <a style="color: #555;" href="https://maayanlab.github.io/appyter/">Appyter Documentation</a><br />
      </p>
    </div>
    <div class="col-md-2 col-sm-3 text-center">
      <a href="https://icahn.mssm.edu/research/bioinformatics" target="_blank">
        <img class="rounded" src="{require('./images/icahn_cb.png')}" style="height: 5rem;">
      </a>
    </div>
    <div class="col-md-2 col-sm-3 text-center">
      <a href="http://lincs-dcic.org" target="_blank">
        <img class="rounded" src="{require('./images/dcic_light.png')}" style="height: 5rem;">
      </a>
    </div>
    <div class="col-md-2 col-sm-3 text-center">
      <a href="https://labs.icahn.mssm.edu/maayanlab/" target="_blank">
        <img class="rounded" src="{require('./images/maayanlab_logo.png')}" style="height: 5rem;">
      </a>
    </div>
    <div class="col-md-2 col-sm-3">
      <div class="row my-2">
        <a class="badge badge-secondary px-2" href="https://github.com/MaayanLab/appyter-catalog" target="_blank">
          <span class="badge badge-light">
            <img src="{require('./images/GitHub-Mark.png')}" style="width: 1rem;">
          </span>
          View source code
        </a>
      </div>
      <div class="row my-2">
        <a class="badge badge-secondary px-2" href="https://github.com/MaayanLab/appyter-catalog/issues/new" target="_blank">
          <span class="badge badge-light">
            <img src="{require('./images/GitHub-Mark.png')}" style="width: 1rem;">
          </span>
          Submit an issue
        </a>
      </div>
    </div>
  </div>
</div>
