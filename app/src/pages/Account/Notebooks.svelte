<script>
  import auth from '../../stores/keycloak_auth_store'
  import Loader from '../../fragments/Loader.svelte'

  const base_url = window.location.origin

  let notebooks

  async function load_notebooks() {
    notebooks = undefined
    const res = await fetch(`${base_url}/postgrest/user_instance`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${$auth.keycloak.token}`,
      },
    })
    notebooks = await res.json()
  }

  async function delete_notebook(id) {
    const res = await fetch(`${base_url}/postgrest/user_instance?id=${encodeURIComponent(`eq.${id}`)}`, {
      method: 'DELETE',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${$auth.keycloak.token}`,
      },
    })
    if (res.status === 204) {
      await load_notebooks()
    } else {
      console.log(await res.text())
    }
  }

  if ($auth.state === 'auth') {
    load_notebooks().catch(e => console.error(e))
  }
</script>

<div class="container flex-grow-1">
  <div class="d-flex justify-content-end py-2">
    <button
      class="btn bg-primary text-white"
      on:click={() => load_notebooks().catch(e => console.error(e))}
    >Refresh</button>
  </div>

  <table class="table table-striped table-fixed">
    <tr>
      <th scope="col">Instance ID</th>
      <th scope="col">Appyter</th>
      <th scope="col">Created</th>
      <th scope="col">Actions</th>
    </tr>
    {#if notebooks === undefined}
      <tr>
        <td class="text-center" colspan="100%"><Loader /></td>
      </tr>
    {:else if notebooks.length === 0}
      <tr>
        <td class="text-center" colspan="100%">No notebooks found</td>
      </tr>
    {:else}
      {#each notebooks as notebook}
        <tr>
          <td><a href="/storage/input/{notebook.file}" download={notebook.filename}>{notebook.filename}</a></td>
          <td>{notebook.ts}</td>
          <td class="text-center"><button
            class="btn bg-danger text-white"
            on:click={() => {
              delete_notebook(notebook.id).catch(e => console.error(e))
            }}
          >Delete</button></td>
        </tr>
      {/each}
    {/if}
  </table>

  <div class="alert alert-primary">
    These are notebooks you've saved. Deleting a notebook from this list will unlink it from your account it will then be subject for deletion as long as no other users have it linked.
  </div>
</div>
