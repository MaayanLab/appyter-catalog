<script>
  import auth from '../../stores/keycloak_auth_store'
</script>

<div class="container flex-grow-1">
  <h1>Welcome, {$auth.keycloak.tokenParsed.given_name}!</h1>
  <p class="lead">Your Appyter user account allows you to keep track of your generated notebooks, and configure custom providers for appyter storage and operation. Navigate with the sections on the left to get started.</p>
</div>