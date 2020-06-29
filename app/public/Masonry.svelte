<script>
  import Macy from 'macy'
  import { onMount, afterUpdate } from 'svelte'

  const id = 'macy-' + Math.random().toString(36).substr(2, 9)
  let masonry

  onMount(() => {
    masonry = new Macy({
      container: '#' + id,
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
          940: 3,
          520: 2,
          400: 1
      },
    })
  })

  afterUpdate(() => masonry.reInit())
</script>

<div id={id}>
  <slot></slot>
</div>
