<script>
  import hash from '../stores/url_hash_store';
  let tab
  $: tab = $hash.params.tab || 'home'
</script>

<div class="container content flex-grow">
  <h1>Account</h1>
  <div class="d-flex">
    <div class="nav flex-column nav-pills m-2" id="v-pills-tab" role="tablist" aria-orientation="vertical">
      <a class="nav-link"
        on:click={() => $hash.params.tab = 'home'}
        class:active={tab === 'home'}
        href="javascript:"
        data-toggle="pill"
        role="tab"
        aria-selected="true">Home</a>
      <a class="nav-link"
        on:click={() => $hash.params.tab = 'notebooks'}
        class:active={tab === 'notebooks'}
        href="javascript:"
        data-toggle="pill"
        role="tab"
        aria-selected="true">Notebooks</a>
      <a class="nav-link"
        on:click={() => $hash.params.tab = 'uploads'}
        class:active={tab === 'uploads'}
        href="javascript:"
        data-toggle="pill"
        role="tab"
        aria-selected="true">Uploads</a>
      <a class="nav-link"
        on:click={() => $hash.params.tab = 'config'}
        class:active={tab === 'config'}
        href="javascript:"
        data-toggle="pill"
        role="tab"
        aria-selected="true">Configuration</a>
    </div>
    <div class="tab-content flex-grow">
      <div class="tab-pane fade show"
        on:click={() => $hash.params.tab = 'home'}
        class:active={tab === 'home'}
        role="tabpanel">
        <p>Welcome to your Appyter user account. This allows you to keep track of your generated notebooks, and configure custom providers for appyter storage and operation.</p>
      </div>
      <div class="tab-pane fade show"
        on:click={() => $hash.params.tab = 'notebooks'}
        class:active={tab === 'notebooks'}
        role="tabpanel">
        <table class="table table-striped">
          <tr>
            <th scope="col">Instance ID</th>
            <th scope="col">Appyter</th>
            <th scope="col">Created</th>
            <th scope="col">Actions</th>
          </tr>
          <tr>
            <td scope="row"><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">6cc2a2c87...</a></td>
            <td><a href="#/Independent_Enrichment_Analysis">Independent Enrichment Analysis</a></td>
            <td>2022-02-24 10:36pm</td>
            <td><a href="javascript:">Delete</a></td>
          </tr>
          <tr>
            <td scope="row"><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">947dbb7c0...</a></td>
            <td><a href="#/Independent_Enrichment_Analysis">Independent Enrichment Analysis</a></td>
            <td>2022-02-24 10:48pm</td>
            <td><a href="javascript:">Delete</a></td>
          </tr>
        </table>
        <p>These are notebooks you've saved. Deleting a notebook from this list will only delete it if no other user has saved it.</p>
      </div>
      <div class="tab-pane fade show"
        on:click={() => $hash.params.tab = 'uploads'}
        class:active={tab === 'uploads'}
        role="tabpanel">
        <table class="table table-striped">
          <tr>
            <th scope="col">Upload ID</th>
            <th scope="col">Filename</th>
            <th scope="col">Created</th>
            <th scope="col">Actions</th>
          </tr>
          <tr>
            <td scope="row"><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">d98b7d978...</a></td>
            <td><a href="">DrugSet_L1000FWD_Signature_Up.txt</a></td>
            <td>2022-02-24 10:30pm</td>
            <td><a href="javascript:">Delete</a></td>
          </tr>
          <tr>
            <td scope="row"><a href="/Independent_Enrichment_Analysis/6cc2a2c87d98b7d978521da5e67947dbb7c09ea1">521da5e67...</a></td>
            <td><a href="">DrugSet_L1000FWD_Signature_Down.txt</a></td>
            <td>2022-02-24 10:31pm</td>
            <td><a href="javascript:">Delete</a></td>
          </tr>
        </table>
        <p>These are data files you've uploaded. Deleting a data file from this list will only delete it if no other user also has it saved.</p>
      </div>
      <div class="tab-pane fade show"
        on:click={() => $hash.params.tab = 'config'}
        class:active={tab === 'config'}
        role="tabpanel">
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
