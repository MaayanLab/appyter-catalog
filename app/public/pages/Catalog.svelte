<script>
  import { onMount } from 'svelte'
  import * as JsSearch from 'js-search'
  import mdIt from 'markdown-it'

  import Home from '../fragments/Home.svelte'
  import Landing from '../fragments/Landing.svelte'
  import SearchBox from '../fragments/SearchBox.svelte'

  import { hash } from '../stores'
  import { hashCode, intToRGB, set_gte } from '../utils.js'

  const base_url = window.location.origin

  // store appyters as list and lookup table based on name slugs
  let { appyters: appyterList } = require('../appyters.json')
  // assemble appyter lookup table
  let appyterLookup = {}
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
    // save a reference in the lookup table
    appyterLookup[name] = appyter
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
    appyterList.sort((a, b) => (b.runs||0) - (a.runs||0))
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
  let loaded = false
  onMount(async () => {
    try {
      await get_pagehits()
    } catch (e) {
      console.error(e)
    }
    loaded = true
  })
</script>

{#if appyter === undefined}
  <div class="row m-0 p-0">
    <div class="col-sm-12 p-0">
      &nbsp;
    </div>
    {#if search}
      <SearchBox />
    {/if}
  </div>
{/if}

<div class="container content flex-grow">
  {#if appyter === undefined}
    <Home
      base_url={base_url}
      loaded={loaded}
      searchAppyters={searchAppyters}
      appyterList={appyterList}
      search={search}
    />
  {:else}
    <Landing
      base_url={base_url}
      appyter={appyter}
    />
  {/if}
</div>
