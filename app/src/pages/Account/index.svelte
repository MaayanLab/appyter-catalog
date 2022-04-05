<script>
  import hash from '@/stores/url_hash_store'
  import auth from '@/stores/keycloak_auth_store'
  import Loader from '@/fragments/Loader.svelte'
  import ErrorOccurred from '@/pages/ErrorOccurred.svelte'

  const pages = {
    '/': {
      name: 'Home',
      component: () => import('./Home.svelte'),
    },
    '/notebooks': {
      name: 'Notebooks',
      component: () => import('./Notebooks.svelte'),
    },
    '/uploads': {
      name: 'Uploads',
      component: () => import('./Uploads.svelte'),
    },
    '/integrations': {
      name: 'Integrations',
      component: () => import('./Integrations/index.svelte'),
    },
    '/account': {
      name: 'Account Management',
      onclick: () => $auth.keycloak.accountManagement(),
    },
    '/logout': {
      name: 'Logout',
      onclick: () => $auth.keycloak.logout(),
    },
    '/404': {
      component: () => import('../NotFound.svelte'),
    },
    '/500': {
      component: async () => ({ default: ErrorOccurred }),
    },
  }

  let path
  $: if ($hash.path !== undefined) {
    let _path = `/${$hash.path.split('/').slice(2, 3).join('/')}`
    if (!(_path in pages)) _path = '/404'
    if (path !== _path) path = _path
  }

  let page
  $: if (path !== undefined) {
    page = undefined
    pages[path].component()
      .then((mod) => page = mod.default)
      .catch((e) => {
        pages['/500'].component()
          .then((mod) => {
            page = mod.default
          })
      })
  }
</script>

{#if $auth.state === 'init'}
  <div class="container d-flex flex-column">
    <h1>Account</h1>
    <Loader />
  </div>
{:else if $auth.state === 'error'}
  <div class="alert alert-warning">
    <p>Account feature is currently unavailable, please try again later.</p>
  </div>
{:else if $auth.state !== 'auth'}
  <div class="alert alert-primary">
    Please login to use the account feature.
  </div>
  <button
    type="button"
    class="btn btn-success btn-lg"
    on:click={() => {
      $auth.keycloak.login()
    }}
  >Login</button>
  <button
    type="button"
    class="btn btn-primary btn-lg"
    on:click={() => {
      $auth.keycloak.register()
    }}
  >Register</button>
{:else}
  <div class="container">
  <div class="d-flex align-items-start">
    <div class="nav flex-column nav-pills me-3" role="tablist" aria-orientation="vertical">
      {#each Object.keys(pages) as p}
        {#if pages[p].name !== undefined}
          <button class="nav-link"
            on:click={() => {
              if (pages[p].onclick !== undefined) {
                pages[p].onclick()
              } else {
                $hash.path = `/account${p}`
              }
            }}
            class:active={path === p}
            data-toggle="pill"
            role="tab"
          >{pages[p].name}</button>
        {/if}
      {/each}
    </div>
    <div class="tab-content">
      <div class="tab-pane show active"
        role="tabpanel">
        {#if page !== undefined}
          <svelte:component this={page} />
        {:else}
          <Loader />
        {/if}
      </div>
    </div>
  </div>
  </div>
{/if}
