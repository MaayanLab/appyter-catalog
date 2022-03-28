import appyterJson from '../../public/appyters.json'
import { writable } from 'svelte/store'
import hash from './url_hash_store.js'

function keycloak_auth_store(keycloakProps) {
  const initStore = {}
  try {
    Object.assign(initStore, {
      state: 'init',
      keycloak: new Keycloak(keycloakProps),
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
        redirectUri: window.location.href + (window.location.href.includes('?') ? '' : '?'),
      }).then(authenticated => {
        keycloak.onTokenExpired = () => {
          console.debug('refreshing expired token...')
          keycloak.updateToken()
            .success(() => {
              set({ state: 'auth', keycloak })
            })
            .error(e => {
              set({ state: 'error', keycloak })
            })
        }
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
        // cleanup keycloak auth params
        hash.update($hash => {
          const params = { ...$hash.params }
          if ('code' in params) delete params['code']
          if ('session_state' in params) delete params['session_state']
          if ('state' in params) delete params['state']
          return { ...$hash, params }
        })
      }).catch(err => {
        console.error(err)
        set({ state: 'error', keycloak })
      })
    }
  })()
  return { subscribe }
}

const auth = keycloak_auth_store(appyterJson.keycloak)

export default auth