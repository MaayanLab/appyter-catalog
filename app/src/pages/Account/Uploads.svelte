<script>
  import auth from '@/stores/keycloak_auth_store'
  import Loader from '@/fragments/Loader.svelte'
  import url_for from '@/utils/url_for'
  import human_size from '@/utils/human_size'

  const base_url = window.location.origin

  let uploads
  let offset = 0
  let limit = 10
  let count

  async function load_uploads({ offset: _offset, limit: _limit }) {
    uploads = undefined
    const res = await fetch(url_for({
      path: `${base_url}/postgrest/user_file`,
      params: {
        offset: _offset,
        limit: _limit
      },
    }), {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Prefer': 'count=estimated',
        'Authorization': `Bearer ${await $auth.keycloak.getValidToken()}`,
      },
    })
    offset = _offset
    limit = _limit
    uploads = await res.json()
    if (uploads.length < limit) {
      count = offset + uploads.length
    } else {
      count = Number(res.headers.get('Content-Range').split('/')[1])
    }
  }

  async function delete_upload(id) {
    const res = await fetch(url_for({
      path: `${base_url}/postgrest/user_file`,
      params: { id: `eq.${id}` }
    }), {
      method: 'DELETE',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${await $auth.keycloak.getValidToken()}`,
      },
    })
    if (res.status === 204) {
      if (offset >= count-1) {
        // deletion would cause this page to disappear, move back a page
        offset = Math.max(0, offset - limit)
      }
      await load_uploads({ offset, limit })
    } else {
      console.error(await res.text())
    }
  }

  if ($auth.state === 'auth') {
    load_uploads({ offset, limit }).catch(e => console.error(e))
  }
</script>

<div class="container flex-grow-1">
  <h1>Uploads</h1>

  <div class="d-flex justify-content-end py-2">
    <button
      class="btn bg-primary text-white"
      on:click={() => load_uploads({ offset, limit }).catch(e => console.error(e))}
    >Refresh</button>
  </div>

  <table class="table table-striped">
    <tr>
      <th scope="col" class="text-left">Filename</th>
      <th scope="col" class="text-left">Size</th>
      <th scope="col" class="text-left">Created</th>
      <th scope="col" class="text-center">Actions</th>
    </tr>
    {#if uploads === undefined}
      <tr>
        <td class="text-center" colspan="100%"><Loader /></td>
      </tr>
    {:else if uploads.length === 0}
      <tr>
        <td class="text-center" colspan="100%">No uploads found</td>
      </tr>
    {:else}
      {#each uploads as upload}
        <tr>
          <td>
            {#if upload.file.id.startsWith('storage://')}
              <a href={upload.file.id.replace(/^storage:\/\//, '/storage/appyters/')} download={upload.filename}>{upload.filename}</a>
            {:else}
              {upload.file.id}
            {/if}
          </td>
          <td>{human_size(upload.file.metadata.size)}</td>
          <td>{upload.ts}</td>
          <td class="text-center"><button
            class="btn btn-sm bg-danger text-white"
            on:click={() => {
              delete_upload(upload.id).catch(e => console.error(e))
            }}
          >Delete</button></td>
        </tr>
      {/each}
      <tr>
        <td colspan="100%">
          {#if !count}
            <div class="text-center">
              No uploads found.
            </div>
          {:else}
            <ul class="pagination">
              <li class="page-item">
                <button
                  class="btn page-link" style="background-color: inherit;"
                  aria-label="Previous"
                  class:disabled={offset === 0}
                  on:click={evt => {
                    load_uploads({ offset: Math.max(0, offset - limit), limit }).catch(e => console.error(e))
                  }}
                >
                  <span aria-hidden="true">&laquo;</span>
                </button>
              </li>
              <li class="page-item">
                <span class="page-link text-black" style="background-color: inherit; border: 0;" aria-label="Next">
                  Showing uploads {offset+1} - {offset + uploads.length} {#if count}of {count}{/if}
                </span>
              </li>
              <li class="page-item">
                <button
                  class="btn page-link" style="background-color: inherit;"
                  class:disabled={offset + limit >= count}
                  on:click={evt => {
                    load_uploads({ offset: Math.min(count - 1, offset + limit), limit }).catch(e => console.error(e))
                  }}
                  aria-label="Next"
                >
                  <span aria-hidden="true">&raquo;</span>
                </button>
              </li>
            </ul>
          {/if}
        </td>
      </tr>
    {/if}
  </table>

  <div class="alert alert-primary">
    These are data files you've uploaded. Deleting a data file from this list will unlink it from your account it will then be subject for deletion as long as no other users have it linked.
  </div>
</div>
