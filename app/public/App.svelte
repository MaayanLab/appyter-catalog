<script>
  import * as JsSearch from 'js-search'

  import mdIt from 'markdown-it'
  const md = mdIt()

  // store templates as list and lookup table based on name slugs
  import templateList from './templates.json'
  const templateLookup = {}
  for (const template of templateList) {
    const {name, description, long_description, ..._} = template
    // modify templates in-place
    Object.assign(template, {
      description: md.render(description || ''),
      long_description: md.render((long_description || description).split('\n').slice(1).join('\n')),
    })
    // save a reference in the lookup table
    templateLookup[name] = template
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
  search.addDocuments(templateList)

  // facilitate search
  let searchString = ''
  const searchTemplates = (searchString) => {
    if (searchString === '') {
      return templateList
    } else {
      return search.search(searchString)
    }
  }

  // sync template variable and url hash
  let template
  const updatehash = () => {
    template = templateLookup[`${window.location.hash || '#'}`.slice(1)]
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
    <h1>Jupyter Template Catalog</h1>
    <h3 class="card-subtitle mb-2 text-muted">A catalog of jupyter templates</h3>
    <hr />
  </div>
</div>
<div class="container">
  {#if template === undefined}
    <div class="row">
      <div class="offset-sm-2 col-sm-8 text-center">
        <input
          type="text"
          class="form-control"
          placeholder="Search templates..."
          aria-label="Search templates"
          bind:value={searchString}
        />
        <p>&nbsp;</p>
      </div>
    </div>
    <div class="row">
      {#each searchTemplates(searchString) as template}
        <div class="col-sm-12 col-md-6 col-xl-4">
          <div class="card">
            <div class="card-body">
              <h5 class="card-title">{template.title}</h5>
              <h6 class="card-subtitle mb-2 text-muted">
                <span class="badge badge-success">v{template.version}</span>
                <span class="badge badge-secondary">{template.license}</span>
                {#each template.tags as tag}
                  <span class="badge badge-primary">{tag}</span>
                {/each}
              </h6>
              <p class="card-text">{@html template.description}</p>
              <a href="#{template.name}" class="btn btn-primary btn-sm stretched-link">
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
        <h3>{template.title}</h3>
        <span class="card-subtitle mb-2 text-muted">
          <span class="badge badge-success">v{template.version}</span>
          <span class="badge badge-secondary">{template.license}</span>
          {#each template.tags as tag}
            <span class="badge badge-primary">{tag}</span>
          {/each}
        </span>
        {#if template.url !== undefined}
          <p><a href="{template.url}">{template.url}</a></p>
        {/if}
        <p>
          <b>Authors:</b><br />
          {#each template.authors as author}
            <span>{author.name} &lt;<a href="mailto:{author.email}">{author.email}</a>&gt;</span><br />
          {/each}
        </p>
        {@html template.long_description}
        <p>&nbsp;</p>
        <a href="./{template.name}/" class="btn btn-primary">Start Template</a>
      </div>
    </div>
  {/if}
</div>
