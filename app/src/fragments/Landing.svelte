<script>
  import AppyterImage from '@/components/AppyterImage.svelte'

  export let base_url
  export let appyter
</script>

<div class="row">
  <div class="col-sm-12 col-md-6 col-lg-4 align-self-center">
    <AppyterImage base_url={base_url} appyter={appyter} />
  </div>
  <div class="col-sm-12 col-md-6 col-lg-8">
    <h2>{appyter.title}</h2>
    <span class="card-subtitle mb-2 text-muted d-flex flex-row flex-wrap">
      <span class="badge bg-success m-1">v{appyter.version}</span>
      <span class="badge bg-secondary m-1">{appyter.license}</span>
      {#each appyter.tags as tag}
        <span class="badge bg-primary m-1">{tag}</span>
      {/each}
    </span>
    <p>
      <a href="{base_url}/{appyter.name}/" class="btn btn-lg btn-danger mb-1">Start Appyter</a>
      <a href="#/running-appyters/?slug={appyter.name}&appyter_version={appyter.version}&run=webform" class="btn btn-lg btn-secondary mb-1">Run Appyter Locally</a>
    </p>
    <p>
      <b>First Published:</b> {(new Date(appyter.creation_timestamp)).toDateString()}<br />
      <b>Last Updated:</b> {(new Date(appyter.update_timestamp)).toDateString()}<br />
    </p>
    <p>
      <b>Authors:</b><br />
      {#each appyter.authors as author}
        <span>
          {author.name}
          <!-- &lt;<a href="mailto:{author.email}">{author.email}</a>&gt; -->
          <br />
        </span>
      {/each}
    </p>
    {#if appyter.url !== undefined}
      <p><a href="{appyter.url}">{appyter.url}</a></p>
    {/if}
  </div>
  <div class="col-sm-12">
    {@html appyter.long_description_html}
  </div>
  <div class="col-sm-12">
    <a href="{base_url}/{appyter.name}/" class="btn btn-lg btn-danger">Start Appyter</a>
    <a href="#/running-appyters/?slug={appyter.name}&appyter_version={appyter.version}&run=webform" class="btn btn-lg btn-secondary">Run Appyter Locally</a>
  </div>
</div>