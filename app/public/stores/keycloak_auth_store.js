import { readable } from 'svelte/store'

function keycloak_auth_store({ url, realm, clientId }) {
  let keycloak, state
  try {
    keycloak = new Keycloak({ url, realm, clientId })
    state = 'init'
  } catch (e) {
    console.error(e)
    keycloak = {}
    state = 'error'
  }
  const { subscribe } = readable({
    keycloak, state,
  }, async (set) => {
    if ('init' in keycloak) {
      const authenticated = await keycloak.init({
        onLoad: 'check-sso',
        silentCheckSsoRedirectUri: `${window.location.origin}/silent-check-sso.html`,
      })
      state = authenticated ? 'auth' : 'guest'
      set({
        state,
        keycloak: {
          ...keycloak,
          logout: () => {
            keycloak.logout()
            state = 'guest'
          }
        },
      })
    }
  })
  return { subscribe }
}

const auth = keycloak_auth_store({
  url: 'https://keycloak.maayanlab.cloud/auth/',
  realm: 'appyters',
  clientId: 'appyter-catalog',
})

export default auth