<script>
  import mdIt from 'markdown-it'
  const md = mdIt()

  // store appyters as list and lookup table based on name slugs
  import appyterList from './appyters.json'
  // assemble appyter lookup table
  const appyterLookup = {}
  for (const appyter of appyterList) {
    const {name, description, long_description, ..._} = appyter
    // modify appyters in-place
    Object.assign(appyter, {
      description: md.render(description || ''),
      long_description: md.render((long_description || description).split('\n').slice(1).join('\n')),
    })
    // save a reference in the lookup table
    appyterLookup[name] = appyter
  }

  // index documents for search
  import * as JsSearch from 'js-search'
  const search = new JsSearch.Search('name')
  search.addIndex('name')
  search.addIndex('title')
  search.addIndex(['authors', 'name']) // TODO: this doesn't work
  search.addIndex(['authors', 'email']) // TODO: this doesn't work
  search.addIndex('description')
  search.addIndex('long_description')
  search.addIndex('license')
  search.addIndex('tags')
  search.addIndex('url')
  search.addDocuments(appyterList)

  // get appyter hits
  // http://localhost:9000/postgrest/pagehits?select=hits&url=eq.test&limit=1
  const base_url = window.location.origin
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
      } else if(url.startsWith(base_url) + '/') {
        const appyter_name = url.slice(base_url.length + 1)
        if (appyterLookup[appyter_name] !== undefined) {
          Object.assign(appyterLookup[appyter_name], { runs: hits })
        }
      }
    }
    searchString = undefined
    searchString = ''
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

  // facilitate search
  let searchString = ''
  const searchAppyters = (searchString) => {
    if (searchString === undefined) return
    if (searchString === '') {
      return appyterList
    } else {
      return search.search(searchString)
    }
  }

  // sync appyter variable and url hash
  let appyter
  const updatehash = () => {
    appyter = appyterLookup[`${window.location.hash || '#'}`.slice(1)]
    pagehit(appyter)
  }
  window.onhashchange = updatehash
  updatehash()

  // things to do on window load
  window.onload = () => {
    get_pagehits()
  }
</script>

<style>
/* sticky footer */
:global(body) {
  display: flex;
  min-height: 100vh;
  flex-direction: column;
}
.content {
  flex: 1;
}

/* card shadow on hover */
.card {
  box-shadow: 0 0px 0px 0 rgba(0,0,0,0.2);
  transition: 0.2s;
}
.card:hover {
  box-shadow: 0 6px 16px 0 rgba(0,0,0,0.2);
}
</style>

<div class="row" style="margin: 0">
  <div class="col-sm-12 text-center">
    <h1>Appyters</h1>
    <h3 class="card-subtitle mb-2 text-muted">A catalog of appyter notebooks</h3>
    <hr />
  </div>
</div>
<div class="container content">
  {#if appyter === undefined}
    <div class="row">
      <div class="offset-sm-2 col-sm-8 text-center">
        <input
          type="text"
          class="form-control"
          placeholder="Search appyters..."
          aria-label="Search appyters"
          bind:value={searchString}
        />
        <p>&nbsp;</p>
      </div>
    </div>
    <div class="row">
      {#each searchAppyters(searchString) as appyter}
        <div class="col-sm-12 col-md-6 col-xl-4">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">{appyter.title}</h5>
              <h6 class="card-subtitle mb-2 text-muted">
                <span class="badge badge-success">v{appyter.version}</span>
                &nbsp;<span class="badge badge-secondary">{appyter.license}</span>
                {#each appyter.tags as tag}
                  &nbsp;<span class="badge badge-primary">{tag}</span>
                {/each}
              </h6>
              <h6 style="color: grey">
                {#if appyter.views }
                  Views: {appyter.views}
                {/if}
                {#if appyter.views && appyter.runs }
                  &nbsp;
                {/if}
                {#if appyter.runs }
                  Runs: n
                {/if}
              </h6>
              <p class="card-text">{@html appyter.description}</p>
              <a href="#{appyter.name}" class="btn btn-primary btn-sm stretched-link">
                Select
              </a>
            </div>
          </div>
        </div>
      {/each}
    </div>
  {:else}
    <div class="row">
      <div class="col-sm-12">
        <h3>{appyter.title}</h3>
        <span class="card-subtitle mb-2 text-muted">
          <span class="badge badge-success">v{appyter.version}</span>
          <span class="badge badge-secondary">{appyter.license}</span>
          {#each appyter.tags as tag}
            <span class="badge badge-primary">{tag}</span>
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
        {@html appyter.long_description}
        <p>&nbsp;</p>
        <a href="./{appyter.name}/" class="btn btn-primary">Start Appyter</a>
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