<script>
  import * as JsSearch from 'js-search'

  import mdIt from 'markdown-it'
  const md = mdIt()

  // store appyters as list and lookup table based on name slugs
  import appyterList from './appyters.json'
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

  // facilitate search
  let searchString = ''
  const searchAppyters = (searchString) => {
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
  }
  window.onhashchange = updatehash
  updatehash()
</script>

<style>
/* card shadow on hover */
.card {
  box-shadow: 0 0px 0px 0 rgba(0,0,0,0.2);
  transition: 0.2s;
}
.card:hover {
  box-shadow: 0 6px 16px 0 rgba(0,0,0,0.2);
}
</style>

<div class="row">
  <div class="col-sm-12 text-center">
    <h1>Appyters</h1>
    <h3 class="card-subtitle mb-2 text-muted">A catalog of appyter notebooks</h3>
    <hr />
  </div>
</div>
<div class="container">
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
                <span class="badge badge-secondary">{appyter.license}</span>
                {#each appyter.tags as tag}
                  <span class="badge badge-primary">{tag}</span>
                {/each}
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
