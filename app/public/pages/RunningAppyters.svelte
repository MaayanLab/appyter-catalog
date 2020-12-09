<script>
  import { hash } from '../stores';

  import StaticMarkdown from '../components/StaticMarkdown.svelte'
  import CodeSnippet from '../components/CodeSnippet.svelte'

  // constants
  let docker_org = 'maayanlab'
  let public_url = window.location.origin
  
  // derived values
  let run
  $: {
    if ($hash.params.run === 'webform') {
      run = 'webform'
    } else if ($hash.params.run === 'cli') {
      run = 'cli'
    } else {
      run = 'notebook'
    }
  }

  let slug_validity
  let version_validity
  let docker_tag
  $: {
    let normalized_slug
    if ($hash.params.slug) {
      const m = /^[0-9A-Za-z_-]+$/.exec($hash.params.slug)
      if (m === null) {
        slug_validity = false
      } else {
        slug_validity = true
      }
      normalized_slug = $hash.params.slug.toLowerCase().replaceAll('-', '_')
    } else {
      slug_validity = undefined
    }
    let docker_image
    if (slug_validity === true) {
      docker_image = `appyter-${normalized_slug}`
    } else {
      docker_image = 'appyter'
    }
    if ($hash.params.version) {
      const m = /^(latest|v\d+\.\d+\.\d+)$/.exec($hash.params.version)
      if (m === null) {
        version_validity = false
      } else {
        version_validity = true
      }
    } else {
      version_validity = undefined
    }
    let docker_version
    if (version_validity === true) {
      docker_version = `:${$hash.params.version}`
    } else {
      docker_version = `:latest`
    }
    docker_tag = `${docker_org}/${docker_image}${docker_version}`
  }

  let id_validity
  $: {
    if ($hash.params.id) {
      const m = /^[0-9A-Za-z-]+$/.exec($hash.params.id)
      if (m === null) {
        id_validity = false
      } else {
        id_validity = true
      }
    } else {
      id_validity = undefined
    }
  }
</script>

<div class="container content flex-grow">
  <StaticMarkdown data={require('./RunningAppyters.md')} />
  <div>
    <form>
      <div class="form-row">
        <div class="form-group col-md-6">
          <label for="appyterSlug">Appyter slug</label>
          <input
            bind:value={$hash.params.slug}
            class:is-valid={slug_validity === true}
            class:is-invalid={slug_validity === false}
            type="text"
            class="form-control"
            id="appyterSlug"
            aria-describedby="appyterSlugHelp"
            placeholder="Example: Drugmonizome_ML"
          >
          {#if version_validity === true}
            <div class="valid-feedback">
              Looks good!
            </div>
          {:else}
            <small id="appyterSlugHelp" class="form-text text-muted">As shown in the case-insensitive url-path of an appyter form or landing page.</small>
            <div class="invalid-feedback">
              This doesn't look like a slug..
            </div>
          {/if}
        </div>
        <div class="form-group col-md-6">
          <label for="appyterVersion">Appyter version</label>
          <input
            bind:value={$hash.params.version}
            class:is-valid={version_validity === true}
            class:is-invalid={version_validity === false}
            type="text"
            class="form-control"
            id="appyterVersion"
            aria-describedby="appyterVersionHelp"
            placeholder="Example: v0.0.1"
          >
          {#if version_validity === true}
            <div class="valid-feedback">
              Looks good!
            </div>
          {:else}
            <small id="appyterVersionHelp" class="form-text text-muted">The version of the appyter you wish to use (defaults to latest)</small>
            <div class="invalid-feedback">
              All appyter versions look something like v0.0.1...
            </div>
          {/if}
        </div>
      </div>
      <div class="form-row">
        <div class="form-group col-md-12">
          <label for="appyterId">Output ID</label>
          <input
            bind:value={$hash.params.id}
            class:is-valid={id_validity === true}
            class:is-invalid={id_validity === false}
            type="text"
            class="form-control"
            id="appyterId"
            aria-describedby="appyterIdHelp"
            placeholder="Example: 725c5cde109032eaeef21bb16654fd922e1de130"
          >
          {#if id_validity === true}
            <div class="valid-feedback">
              Looks good!
            </div>
          {:else}
            <small id="appyterIdHelp" class="form-text text-muted">The identifier for your appyter as it appears at the end of the url (it's a large high-entropy alpha-numeric string)</small>
            <div class="invalid-feedback">
              This doesn't look right..
            </div>
          {/if}
        </div>
      </div>
    </form>
  </div>
  <h3>Step 3. Run the appyter</h3>
  <p>The instructions will update in real time as the information provided in Step 2 is updated.</p>
  <ul class="nav nav-tabs">
    <li class="nav-item">
      <a
        on:click={() => $hash.params.run = 'notebook'}
        class:active={run === 'notebook'}
        class="nav-link"
        href="javascript:"
      >Jupyter Notebook</a>
    </li>
    <li class="nav-item">
      <a
        on:click={() => $hash.params.run = 'webform'}
        class:active={run === 'webform'}
        class="nav-link"
        href="javascript:"
      >Web Form</a>
    </li>
    <li class="nav-item">
      <a
        on:click={() => $hash.params.run = 'cli'}
        class:active={run === 'cli'}
        class="nav-link"
        href="javascript:"
      >Command Line Application</a>
    </li>
  </ul>
  <div class="card">
    <div class="card-body">
      <div class="tab-content">
        <div class="tab-pane fade show active" role="tabpanel">
          {#if run === 'notebook'}
            <h4>Running as a jupyter notebook</h4>
            <p>
              Serving the results of an appyter as a live jupyter notebook can be done locally if you install all necessary dependencies,
              which are listed in each appyter's public source code but this can be a challenge. Instead, you can use the appyter's docker container
              which was used to originally execute the notebook.
            </p>
            {#if id_validity === true}
              <p>
                 The following command will fetch your notebook and serve it using a standard jupyter notebook.
              </p>
              <CodeSnippet
                code={`
                  docker run -p 5000:5000 -it ${docker_tag} appyter fetch-and-serve ${public_url}/${$hash.params.slug}/${$hash.params.id}/
                `}
              />
            {:else}
              <p>
                 The following command will serve a jupyter notebook with the same environment of the appyter.
                 If you provide the notebook identifier, we'll update the command to automatically fetch your notebook.
              </p>
              <CodeSnippet
                code={`
                  docker run -p 5000:5000 -it ${docker_tag} appyter serve
                `}
              />
            {/if}
            <p>The notebook should then be accessible in your web browser served from the container at <b>port 5000</b>.</p>
            <p>On Mac OSX and Linux, <a href="http://localhost:5000">http://localhost:5000</a> can typically be used.</p>
            <p>
              On Windows or docker installations which use docker-machine, it may be necessary to use the docker-machine's ip address.
              This can usually be found with the <code>docker-machine ip</code> command.
              For further information, please consult the <a href="https://docs.docker.com/">Docker documentation</a>.
            </p>
          {:else if run === 'webform'}
            <h4>Running as a web-form</h4>
            <p>Executing the appyter as a web-form will allow you to serve the same type of application that we are serving but all execution happens in the docker image.</p>
            <CodeSnippet
              code={`
                docker run -p 5000:5000 -it ${docker_tag}
              `}
            />
            <p>The notebook should then be accessible in your web browser served from the container at <b>port 5000</b>.</p>
            <p>On Mac OSX and Linux, <a href="http://localhost:5000">http://localhost:5000</a> can typically be used.</p>
            <p>
              On Windows or docker installations which use docker-machine, it may be necessary to use the docker-machine's ip address.
              This can usually be found with the <code>docker-machine ip</code> command.
              For further information, please consult the <a href="https://docs.docker.com/">Docker documentation</a>.
            </p>
          {:else if run === 'cli'}
            <h4>Running as a command line application</h4>
            <p>Interacting with appyter as a command line application will allow you to use it as part of a workflow.</p>
            <CodeSnippet
              code={`
                # Inspect an appyter for its available fields
                docker run -it ${docker_tag} appyter nbinspect
                
                # Inspect construct an appyter using your own input data
                docker run -it ${docker_tag} appyter nbconstruct < input.json > output.ipynb
                
                # Execute a constructed jupyter notebook on the CLI
                # exit code 0 for success, 1 for failure
                docker run -it ${docker_tag} appyter nbexecute output.ipynb
                
                # Review the appyter CLI options
                docker run -it ${docker_tag} appyter --help
              `}
            />
            <p>
              Alternatively, you can spawn a bash terminal in the docker container and run the commands directly in the container allowing you to ommit the docker run part.
            </p>
            <CodeSnippet
              code={`
                docker run -it ${docker_tag} /bin/bash
              `}
            />
          {/if}
        </div>
      </div>
    </div>
  </div>
  <p>&nbsp;</p>
</div>
