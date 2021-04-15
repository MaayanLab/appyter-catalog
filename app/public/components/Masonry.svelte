<script>
  import Macy from 'macy'
  import { tick } from 'svelte'

  export let items

  let masonry
  let container

  $: if (container) {
    masonry = new Macy({
      container,
      trueOrder: true,
      waitForImages: false,
      useOwnImageLoader: false,
      mobileFirst: true,
      useContainerForBreakpoints: true,
      margin: {
        y: 16,
        x: '2%',
      },
      breakAt: {
        1000: 3,
        800: 2,
        400: 1
      },
    })
  }

  $: if (masonry && items) {
    tick()
      .then(_ => masonry.recalculate(true, true))
      .then(_ => tick())
      .then(_ => masonry.recalculate(true, true))
  }
</script>

<div bind:this={container}>
  {#each items as item}
    <slot item={item}></slot>
  {/each}
</div>
