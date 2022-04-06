<script>
  export let base_url
  export let appyter

  function localize_appyter_image(base_url, { name, image }) {
    if (/^https?:\/\//.exec(image) !== null) {
      return image
    } else {
      return `${base_url}/${name}/static/${image}`
    }
  }
</script>

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