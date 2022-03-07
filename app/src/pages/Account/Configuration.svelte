<script>
  import auth from '../../stores/keycloak_auth_store'

  const base_url = window.location.origin

  let savedConfig = {
    storage: '',
    executor: '',
  }
  let config = {...savedConfig}

  async function ensure_config(props) {
    const res = await fetch(`${base_url}/postgrest/rpc/user_config`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${$auth.keycloak.token}`,
      },
      body: JSON.stringify(props !== undefined ? props : { config: null }),
    })
    if (res.status === 200) {
      savedConfig = await res.json()
      config = {...savedConfig}
    } else {
      console.error(await res.text())
    }
    // TODO: handle errors
  }

  $: if ($auth.state === 'auth') {
    ensure_config().catch(e => console.error(e))
  }
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
      bind:value={config.storage}
    >
  </div>
  <div class="form-group col-sm-12">
    <label for="executor">Custom Execution Protocol</label>
    <input
      type="text"
      class="form-control"
      id="executor"
      placeholder="e.g. wes://my-wes-endpoint"
      bind:value={config.executor}
    >
  </div>
  <button
    type="button"
    class="btn btn-success btn-lg"
    on:click={() => {
      ensure_config({ config })
    }}
    class:disabled={!(config.storage !== savedConfig.storage || config.executor !== savedConfig.executor)}
  >Save</button>
  <button
    type="button"
    class="btn btn-danger btn-lg"
    class:disabled={savedConfig.storage === '' && savedConfig.executor === ''}
    on:click={() => {
      ensure_config({ config: { storage: '', executor: '' } })
    }}
  >Delete</button>
</div>
