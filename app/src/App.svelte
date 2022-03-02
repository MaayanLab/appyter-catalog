<script>
  import hash from './stores/url_hash_store'

  import Header from './fragments/Header.svelte'
  import Footer from './fragments/Footer.svelte'
  import Loader from './fragments/Loader.svelte'
  import ErrorOccurred from './pages/ErrorOccurred.svelte'

  let pages = {
    '': () => import('./pages/Catalog.svelte'),
    '/what-is-an-appyter': () => import('./pages/WhatIsAnAppyter.svelte'),
    '/creating-appyters': () => import('./pages/CreatingAppyters.svelte'),
    '/publishing-appyters': () => import('./pages/PublishingAppyters.svelte'),
    '/running-appyters': () => import('./pages/RunningAppyters.svelte'),
    '/about': () => import('./pages/About.svelte'),
    '/account': () => import('./pages/Account.svelte'),
    '/404': () => import('./pages/NotFound.svelte'),
    '/500': async () => ({ default: ErrorOccurred }),
  }
  let path
  $: if ($hash.path !== undefined) {
    let _path = $hash.path.split('/').slice(0, -1).join('/')
    if (!(_path in pages)) _path = '/404'
    if (path !== _path) path = _path
  }

  let page
  $: if (path !== undefined) {
    page = undefined
    pages[path]()
      .then((mod) => page = mod.default)
      .catch((e) => {
        pages['/500']()
          .then((mod) => {
            page = mod.default
          })
      })
  }
</script>

<style>
  /* sticky footer */
  :global(body) {
    display: flex;
    min-height: 100vh;
    flex-direction: column;
    min-width: 540px;
  }

  :global(body) {
    background-color: #f5f5f5;
  }

  :global(a) {
    text-decoration: none;
  }
  :global(button), :global(button:hover), :global(button:focus), :global(button:active) {
    border: none !important;
    outline: none !important;
  }
</style>

<Header />
<div class="container flex-grow-1">
  {#if page !== undefined}
    <svelte:component this={page} />
  {:else}
    <Loader />
  {/if}
</div>
<Footer />
