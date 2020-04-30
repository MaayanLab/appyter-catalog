<script>
  // store templates as list and lookup table based on name slugs
  import templateList from './templates.json'
  const templateLookup = {}
  for (const template of templateList) {
    templateLookup[template.name] = template
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
      {#each templateList as template}
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
              <p class="card-text">{template.description}</p>
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
        <p>&nbsp;</p>
        <p>{template.description}</p>
        {#if template.url !== undefined}
          <a href="{template.url}">{template.url}</a>
        {/if}
        <p>&nbsp;</p>
        <p>
          <b>Authors:</b><br />
          {#each template.authors as author}
            <span>{author.name} &lt;<a href="mailto:{author.email}">{author.email}</a>&gt;</span><br />
          {/each}
        </p>
        <a href="./{template.name}/" class="btn btn-primary">Start Template</a>
      </div>
    </div>
  {/if}
</div>
