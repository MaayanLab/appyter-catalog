<script>
  import auth from '@/stores/keycloak_auth_store'
  import cavaticaGuideApiKey from '@/public/images/CAVATICA-guide-apikey.png'
  import cavaticaGuideProject from '@/public/images/CAVATICA-guide-project.png'

  const base_url = window.location.origin
  let savedConfig = {
    cavatica_api_key: '',
    cavatica_project: '',
  }
  let config = {}

  async function ensure_config(props) {
    const res = await fetch(`${base_url}/postgrest/rpc/user_config`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Authorization': `Bearer ${await $auth.keycloak.getValidToken()}`,
      },
      body: JSON.stringify(props !== undefined ? props : { config: null }),
    })
    if (res.status === 200) {
      savedConfig = await res.json()
      config = {...savedConfig}
    } else {
      console.error(await res.text())
    }
  }

  $: console.log(config)
  $: if ($auth.state === 'auth') {
    ensure_config().catch(e => console.error(e))
  }
</script>

<div class="container flex-grow-1">
  <h1>CAVATICA Integration</h1>

  <p><a href="https://www.cavatica.org/">CAVATICA</a> is a data analysis and sharing platform designed to accelerate discovery in a scalable, cloud-based compute environment where data, results, and workflows are shared among the world's research community. Developed by Seven Bridges and funded in-part by a grant from the National Institutes of Health (NIH) Common Fund, CAVATICA is continuously updated with new tools and datasets.</p>

  <p>CAVATICA offers secure storage and compute in a cloud environment. Appyter integration with CAVATICA enables you to execute Appyters against data in a CAVATICA project and use CAVATICA-managed computational resources.</p>

  <p>To use CAVATICA, you must <a href="https://pgc-accounts.sbgenomics.com/auth/login" target="_blank">login or create an account</a>, then a project should be established which will be used for file storage and executions that are configured.</p>
  <img
    class="img-fluid"
    src={cavaticaGuideProject}
    alt="CAVATICA Project Creation Guide" />

  <p>To use CAVATICA you must register your CAVATICA API Key with Appyters, this key can be located at <a href="https://cavatica.sbgenomics.com/developer/token" target="_blank">https://cavatica.sbgenomics.com/developer/token</a></p>
  <img
    class="img-fluid"
    src={cavaticaGuideApiKey}
    alt="CAVATICA API Key Guide" />

  <div class="container flex-grow-1 my-4">
    <div class="form-group col-sm-12">
      <label for="cavatica_api_key">CAVATICA API Key</label>
      <input
        type="text"
        class="form-control"
        id="cavatica_api_key"
        placeholder="e.g. 08cd35123..."
        bind:value={config.cavatica_api_key}
      />
    </div>
    <div class="form-group col-sm-12">
      <label for="cavatica_project">CAVATICA Project (default)</label>
      <input
        type="text"
        class="form-control"
        id="cavatica_project"
        placeholder="e.g. youruser/yourproject"
        bind:value={config.cavatica_project}
      />
    </div>
  </div>
    <button
      type="button"
      class="btn btn-success btn-lg"
      on:click={() => {
        ensure_config({ config })
      }}
      class:disabled={!(JSON.stringify(config) !== JSON.stringify(savedConfig))}
    >Save</button>
    <button
      type="button"
      class="btn btn-danger btn-lg"
      class:disabled={savedConfig.cavatica_api_key === '' && savedConfig.cavatica_project === ''}
      on:click={() => {
        ensure_config({ config: { cavatica_api_key: '', cavatica_project: '' } })
      }}
    >Delete</button>

  <p>&nbsp;</p>

  <div class="alert alert-warning">
    CAVATICA Integration support is currently in BETA, we appreciate any feedback you can provide. See <code>Contact Us</code> or <code>Submit an Issue</code> links at the bottom of this page.
  </div>
</div>