<script>
  import hash from '../stores/url_hash_store';
  import auth from '../stores/keycloak_auth_store'
  import Loader from '../fragments/Loader.svelte'
</script>

{#if $auth.state === 'init'}
  <div class="container d-flex flex-column">
    <h1>Account</h1>
    <Loader />
  </div>
{:else if $auth.state === 'error'}
  <div class="alert alert-warning">
    <p>Account feature is currently unavailable, please try again later.</p>
  </div>
{:else if $auth.state !== 'auth'}
  <div class="alert alert-primary">
    Please login to use the account feature.
  </div>
  <button
    type="button"
    class="btn btn-success btn-lg"
    on:click={() => {
      $auth.keycloak.login()
    }}
  >Login</button>
  <button
    type="button"
    class="btn btn-primary btn-lg"
    on:click={() => {
      $auth.keycloak.register()
    }}
  >Register</button>
{:else}
  <div class="container">
  <div class="d-flex align-items-start">
    <div class="nav flex-column nav-pills me-3" role="tablist" aria-orientation="vertical">
      <button class="nav-link"
        on:click={() => $hash.params.tab = undefined}
        class:active={!$hash.params.tab}
        data-toggle="pill"
        role="tab"
      >Home</button>
      <button class="nav-link"
        on:click={() => $hash.params.tab = 'notebooks'}
        class:active={$hash.params.tab === 'notebooks'}
        data-toggle="pill"
        role="tab"
      >Notebooks</button>
      <button class="nav-link"
        on:click={() => $hash.params.tab = 'uploads'}
        class:active={$hash.params.tab === 'uploads'}
        data-toggle="pill"
        role="tab"
      >Uploads</button>
      <button class="nav-link"
        on:click={() => $hash.params.tab = 'config'}
        class:active={$hash.params.tab === 'config'}
        data-toggle="pill"
        role="tab"
      >Configuration</button>
      <button class="nav-link"
        on:click={() => $auth.keycloak.accountManagement()}
        data-toggle="pill"
        role="tab"
        aria-selected="true">Account Management</button>
      <button class="nav-link"
        on:click={() => $auth.keycloak.logout()}
        data-toggle="pill"
        role="tab"
        aria-selected="true">Logout</button>
    </div>
    <div class="tab-content">
      <div class="tab-pane fade"
        on:click={() => $hash.params.tab = undefined}
        class:show={!$hash.params.tab}
        class:active={!$hash.params.tab}
        role="tabpanel">
        <div class="container flex-grow-1">
          <h1 class="display-5">Welcome, {$auth.keycloak.tokenParsed.given_name}!</h1>
          <p class="lead">Your Appyter user account allows you to keep track of your generated notebooks, and configure custom providers for appyter storage and operation. Navigate with the sections on the left to get started.</p>
        </div>
      </div>
      <div class="tab-pane fade"
        on:click={() => $hash.params.tab = 'notebooks'}
        class:show={$hash.params.tab === 'notebooks'}
        class:active={$hash.params.tab === 'notebooks'}
        role="tabpanel">
        <div class="container flex-grow-1">
          <table class="table table-striped">
            <tr>
              <th scope="col">Instance ID</th>
              <th scope="col">Appyter</th>
              <th scope="col">Created</th>
              <th scope="col">Actions</th>
            </tr>
            <tr>
              <td><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">6cc2a2c87...</a></td>
              <td><a href="#/Independent_Enrichment_Analysis">Independent Enrichment Analysis</a></td>
              <td>2022-02-24 10:36pm</td>
              <td><button>Delete</button></td>
            </tr>
            <tr>
              <td><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">947dbb7c0...</a></td>
              <td><a href="#/Independent_Enrichment_Analysis">Independent Enrichment Analysis</a></td>
              <td>2022-02-24 10:48pm</td>
              <td><button>Delete</button></td>
            </tr>
          </table>
          <p>These are notebooks you've saved. Deleting a notebook from this list will only delete it if no other user has saved it.</p>
        </div>
      </div>
      <div class="tab-pane fade"
        on:click={() => $hash.params.tab = 'uploads'}
        class:show={$hash.params.tab === 'uploads'}
        class:active={$hash.params.tab === 'uploads'}
        role="tabpanel">
        <div class="container flex-grow-1">
          <table class="table table-striped">
            <tr>
              <th scope="col">Upload ID</th>
              <th scope="col">Filename</th>
              <th scope="col">Created</th>
              <th scope="col">Actions</th>
            </tr>
            <tr>
              <td><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">d98b7d978...</a></td>
              <td><a href="DrugSet_L1000FWD_Signature_Up.txt">DrugSet_L1000FWD_Signature_Up.txt</a></td>
              <td>2022-02-24 10:30pm</td>
              <td><button>Delete</button></td>
            </tr>
            <tr>
              <td><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">521da5e67...</a></td>
              <td><a href="DrugSet_L1000FWD_Signature_Down.txt">DrugSet_L1000FWD_Signature_Down.txt</a></td>
              <td>2022-02-24 10:31pm</td>
              <td><button>Delete</button></td>
            </tr>
          </table>
          <p>These are data files you've uploaded. Deleting a data file from this list will only delete it if no other user also has it saved.</p>
        </div>
      </div>
      <div class="tab-pane fade"
        on:click={() => $hash.params.tab = 'config'}
        class:show={$hash.params.tab === 'config'}
        class:active={$hash.params.tab === 'config'}
        role="tabpanel">
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
      </div>
    </div>
  </div>
  </div>
{/if}
