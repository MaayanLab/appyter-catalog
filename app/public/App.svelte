<script>
  import { hash } from './stores'

  import Header from './fragments/Header.svelte'
  import Footer from './fragments/Footer.svelte'

  let pages = {
    '': require('./pages/Catalog.svelte').default,
    '/what-is-an-appyter': require('./pages/WhatIsAnAppyter.svelte').default,
    '/using-appyters': require('./pages/UsingAppyters.svelte').default,
    '/creating-appyters': require('./pages/CreatingAppyters.svelte').default,
    '/publishing-appyters': require('./pages/PublishingAppyters.svelte').default,
    '/about': require('./pages/About.svelte').default,
    [undefined]: require('./pages/NotFound.svelte').default,
  }
</script>

<style>
  /* sticky footer */
  :global(body) {
    display: flex;
    min-height: 100vh;
    flex-direction: column;
  }

  :global(.flex-grow) {
    flex: 1 0 auto;
  }

  :global(body) {
    background-color: #f5f5f5;
  }
</style>

<Header />
<svelte:component
  this={
    pages[$hash.path.split('/').slice(0, -1).join('/')]
    || pages[undefined]
  }
/>
<Footer />
