<script>
  import Macy from 'macy'
  import { onMount, afterUpdate } from 'svelte'

  let masonry
  let container

  onMount(() => {
    masonry = new Macy({
      container,
      trueOrder: true,
      waitForImages: false,
      useOwnImageLoader: false,
      debug: true,
      mobileFirst: true,
      columns: 1,
      margin: {
        y: 16,
        x: '2%',
      },
      breakAt: {
        1200: 3,
        992: 2,
        576: 1
      },
    })
  })

  afterUpdate(() => {
    if (masonry !== undefined) {
      masonry.recalculate(true)
    }
  })
</script>

<div bind:this={container}>
  <slot></slot>
</div>
