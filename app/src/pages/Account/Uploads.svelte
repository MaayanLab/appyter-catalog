<script>
  import auth from '../../stores/keycloak_auth_store'
  import Loader from '../../fragments/Loader.svelte'

  const base_url = window.location.origin

  let uploads

  async function load_uploads() {
    uploads = undefined
    const res = await fetch(`${base_url}/postgrest/user_file`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${$auth.keycloak.token}`,
      },
    })
    uploads = await res.json()
  }

  async function delete_upload(id) {
    const res = await fetch(`${base_url}/postgrest/user_file?id=${encodeURIComponent(`eq.${id}`)}`, {
      method: 'DELETE',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${$auth.keycloak.token}`,
      },
    })
    if (res.status === 204) {
      await load_uploads()
    } else {
      console.log(await res.text())
    }
  }

  if ($auth.state === 'auth') {
    load_uploads().catch(e => console.error(e))
  }
</script>

<div class="container flex-grow-1">
  <div class="d-flex justify-content-end py-2">
    <button
      class="btn bg-primary text-white"
      on:click={() => load_uploads().catch(e => console.error(e))}
    >Refresh</button>
  </div>

  <table class="table table-striped table-fixed">
    <tr>
      <th scope="col" class="text-left">Filename</th>
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
          <td><a href="/storage/input/{upload.file}" download={upload.filename}>{upload.filename}</a></td>
          <td>{upload.ts}</td>
          <td class="text-center"><button
            class="btn bg-danger text-white"
            on:click={() => {
              delete_upload(upload.id).catch(e => console.error(e))
            }}
          >Delete</button></td>
        </tr>
      {/each}
    {/if}
  </table>

  <div class="alert alert-primary">
    These are data files you've uploaded. Deleting a data file from this list will unlink it from your account it will then be subject for deletion as long as no other users have it linked.
  </div>
</div>
