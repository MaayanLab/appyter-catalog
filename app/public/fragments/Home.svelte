<script>
  import { hash } from '../stores'
  import Masonry from '../components/Masonry'
  import { localize_appyter_image } from '../utils.js'

  export let base_url
  export let loaded
  export let searchAppyters
  export let appyterList
  export let search
</script>

<style>
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

  .text-placeholder {
    height: 8px;
    margin-bottom: 5px;
  }

  .placeholder {
    display: inline-block;
    background-color: #b0c0c7;
    animation-name: shine;
    animation-duration: 2.4s;
    animation-iteration-count: infinite;
  }

  @keyframes shine {
    0% { 
      opacity: 1;
    }
    50% {
      opacity: 0.35;
    }
    100% {
      opacity: 1;
    }
  }
</style>

<Masonry
  items={!loaded ? [1, 2, 3] : searchAppyters(appyterList, search, $hash.params.q, $hash.params.tags)}
  let:item={appyter}
>
  {#if !loaded}
    <div class="card">
      <div
        class="card-img-top placeholder"
        style={[
          `width: 100%`,
          `background-color: grey;`,
          `padding-top: 56.25%`, // 720 / 1280 = 0.5625, this preserves aspect ratio of div
        ].filter((el) => el !== undefined).join('; ')}
      ></div>
      <div class="card-body">
        <h3 class="card-title">
          <div class="text-placeholder w-100"></div>
        </h3>
        <p class="card-text">
          <div class="text-placeholder placeholder w-100"></div>
          <div class="text-placeholder placeholder w-75"></div>
          <div class="text-placeholder placeholder w-100"></div>
        <div class="pb-4 d-flex flex-row flex-wrap">
          <div class="text-placeholder placeholder"></div>
        </div>
        <button class="btn btn-secondary btn-sm disabled">
          View Details
        </button>
        <button class="btn btn-sm btn-danger disabled">
          Start Appyter
        </button>
        <button class="btn btn-sm btn-secondary disabled">
          Run Appyter Locally
        </button>
      </div>
    </div>
  {:else}
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
              `background-image: url('${localize_appyter_image(base_url, appyter)}')`
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
          <button class="btn btn-primary btn-sm mb-1">
            View Details
          </button>
          <a href="{base_url}/{appyter.name}/" class="btn btn-sm btn-danger mb-1">
            Start Appyter
          </a>
          <a href="#/running-appyters/?slug={appyter.name}&appyter_version={appyter.version}&run=webform" class="btn btn-sm btn-secondary mb-1">
            Run Locally
          </a>
        </div>
      </div>
    </a>
  {/if}
</Masonry>