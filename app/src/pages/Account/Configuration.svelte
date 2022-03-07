<script>
  import hash from '../../stores/url_hash_store'
</script>

<div class="container flex-grow-1">
  <p><b>ALPHA</b> For advanced users, managed configuration coming soon.</p>
  <div class="form-group col-sm-12">
    <label for="storage">Custom Storage Protocol</label>
    <input
      type="text"
      class="form-control"
      id="storage"
      placeholder="e.g. s3://my-data-bucket"
      bind:value={$hash.params.storage}
    >
  </div>
  <div class="form-group col-sm-12">
    <label for="executor">Custom Execution Protocol</label>
    <input
      type="text"
      class="form-control"
      id="executor"
      placeholder="e.g. wes://my-wes-endpoint"
      bind:value={$hash.params.executor}
    >
  </div>
  <button
    type="button"
    class="btn btn-success btn-lg"
    on:click={() => {
      // TODO: apply
    }}
    class:disabled={!($hash.params.storage || $hash.params.executor)}
  >Save</button>
  <button
    type="button"
    class="btn btn-danger btn-lg"
    on:click={() => {
      $hash.params.storage = ''
      $hash.params.executor = ''
      // TODO: apply
    }}
  >Delete</button>
</div>
