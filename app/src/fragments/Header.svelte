<script>
  import hash from '../stores/url_hash_store'
  import appyterLogo from '../../public/images/appyters_logo.svg'
</script>

<style>
  a {
    text-decoration: none;
    color: inherit;
  }
  a:hover {
    text-decoration: underline;
    cursor: pointer;
  }
  a.selected {
    font-weight: bold;
  }
</style>

<div class="container-fluid bg-white pb-2 mb-4">
  <div class="container">
    <div class="row">
      <div class="col-md-12 col-lg-4">
        <h1 class="m-0">
          <a href="#/">
            <img
              src="{appyterLogo}"
              class="img-fluid w-100 p-4"
              alt="Appyters"
            />
          </a>
        </h1>
      </div>
      <div class="col-sm-12 col-lg-8 offset-xl-1 col-xl-7 text-center my-4" style="white-space: nowrap">
        <span class="text-muted">
          <a
            href="#/what-is-an-appyter/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/what-is-an-appyter'}
          >
            What is an Appyter?
          </a>
            |
          <a
            href="#/creating-appyters/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/creating-appyters'}
          >
            Creating Appyters
          </a>
           | 
          <a
            href="#/publishing-appyters/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/publishing-appyters'}
          >
            Publishing Appyters
          </a>
           | 
          <a
            href="#/about/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/about'}
          >
            About
          </a>
           | 
          <a
            href="#/account/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/account'}
          >
            Account
          </a>
        </span>
      </div>
    </div>
  </div>
</div>
