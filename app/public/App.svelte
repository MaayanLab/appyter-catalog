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