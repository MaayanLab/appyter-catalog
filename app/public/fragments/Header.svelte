<script>
  import { hash } from '../stores'
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

<div class="container-fluid bg-white">
  <div class="container">
    <div class="row">
      <div class="col-md-12 col-lg-4 col-xl-5">
        <h1>
          <a href="#/">
            <img
              src="{require('../images/appyters_logo.svg')}"
              class="img-fluid w-100 p-2"
              alt="Appyters"
            />
          </a>
        </h1>
      </div>
      <div class="col-sm-12 col-lg-8 col-xl-7 text-right my-4" style="white-space: nowrap">
        <span class="text-muted">
          <a
            href="#/what-is-an-appyter/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/what-is-an-appyter'}
          >
            What is an Appyter?
          </a>
            |
          <a
            href="#/using-appyters/"
            class:selected={$hash.path.split('/').slice(0, -1).join('/')==='/using-appyters'}
          >
            Using Appyters
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
        </span>
      </div>
    </div>
  </div>
</div>
