<script>
  import { localize_appyter_image } from '../utils.js'

  export let base_url
  export let appyter
</script>

<div class="row">
  <div class="col-sm-12 col-md-6 col-lg-4 align-self-center">
    <div 
      class="card-img-top"
      style={[
        `background-color: ${appyter.color}`,
        appyter.image !== undefined ? (
          `background-image: url('${localize_appyter_image(base_url, appyter)}')`
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
    <a href="#/running-appyters/?slug={appyter.name}&appyter_version={appyter.version}&run=webform" class="btn btn-secondary">Run Appyter Locally</a>
  </div>
</div>