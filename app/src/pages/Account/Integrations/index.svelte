<script>
  import hash from '../../../stores/url_hash_store'
  import Loader from '../../../fragments/Loader.svelte'
  import ErrorOccurred from '../../ErrorOccurred.svelte'

  const pages = {
    '/': {
      component: () => import('./Listing.svelte'),
    },
    '/404': {
      component: () => import('../../NotFound.svelte'),
    },
    '/500': {
      component: async () => ({ default: ErrorOccurred }),
    },
  }

  let path
  $: if ($hash.path !== undefined) {
    let _path = `/${$hash.path.split('/').slice(3, 4).join('/')}`
    console.log(_path)
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

{#if page !== undefined}
  <svelte:component this={page} />
{:else}
  <Loader />
{/if}