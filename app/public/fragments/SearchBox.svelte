<script>
  import hash from '../stores/url_hash_store'

  let { global_tags } = require('../appyters.json')
</script>

<div class="container">
  <div class="col-sm-12 col-xl-10 mb-3">
    <div class="row">
      <div class="col-sm-12">
        <input
          type="text"
          class="form-control"
          placeholder="Search appyters..."
          aria-label="Search appyters"
          bind:value={$hash.params.q}
        />
      </div>
      {#if $hash.params.tags}
        <div class="col-sm-12 mb-2">
          {#each $hash.params.tags.split(';') as tag}
            <a
              href="javascript:"
              class="badge badge-success m-1 p-1 text-white"
              on:click={() => $hash.params.tags = $hash.params.tags.split(';').filter(t => t !== tag).join(';')}
            >{tag}</a>
          {/each}
        </div>
      {/if}
      <div class="col-sm-12">
        {#each global_tags.filter(tag => ($hash.params.tags || '').split(';').indexOf(tag) === -1) as tag}
          <a
            href="javascript:"
            class="badge badge-primary m-1 p-1 text-white"
            on:click={() => $hash.params.tags = [...($hash.params.tags || '').split(';').filter(t => t !== '' && t !== tag), tag].join(';')}
          >{tag}</a>
        {/each}
      </div>
    </div>
  </div>
</div>