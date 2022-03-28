<script>
  import cavaticaLogo from '../../../../public/images/CAVATICA.png'
  import hash from '../../../stores/url_hash_store'
</script>

<div class="container flex-grow-1">
  <h1>Integrations</h1>
  <p class="lead">Click on a platform logo to get started with that platform's Appyter integration.</p>
  <div class="row">
    <div class="col-sm-12 col-md-6">
      <button
        on:click={() => {
          $hash.path = $hash.path + `/cavatica`
        }}
      >
        <img class="img-fluid" src={cavaticaLogo} alt="CAVATICA" />
      </button>
    </div>
  </div>
</div>
