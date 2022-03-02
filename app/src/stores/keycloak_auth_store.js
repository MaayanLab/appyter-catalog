import { readable, writable } from 'svelte/store'

function keycloak_auth_store({ url, realm, clientId }) {
  const initStore = {}
  try {
    Object.assign(initStore, {
      state: 'init',
      keycloak: new Keycloak({ url, realm, clientId }),
    })
  } catch (e) {
    console.error(e)
    Object.assign(initStore, {
      state: 'error',
      keycloak: {},
    })
  }
  const { subscribe, set } = writable(initStore)
  ;(async () => {
    const { keycloak } = initStore
    if ('init' in keycloak) {
      keycloak.init({
        onLoad: 'check-sso',
        silentCheckSsoRedirectUri: `${window.location.origin}/silent-check-sso.html`,
      }).then(authenticated => {
        console.log({ authenticated })
        set({
          state: authenticated ? 'auth' : 'guest',
          keycloak: {
            ...keycloak,
            logout: () => {
              keycloak.logout()
              set({ state: 'guest', keycloak })
            }
          },
        })
      }).catch(err => {
        console.error(err)
        set({ state: 'error', keycloak })
      })
    }
  })()
  return { subscribe }
}

const auth = keycloak_auth_store({
  url: 'https://keycloak.maayanlab.cloud/auth/',
  realm: 'appyters',
  clientId: 'appyter-catalog',
})

export default auth