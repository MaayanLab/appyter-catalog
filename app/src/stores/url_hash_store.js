import { writable } from 'svelte/store'

function params_parse(params_encoded) {
  if (params_encoded === '') return {}
  return params_encoded
    .split('&')
    .map(param_pair => param_pair.split('=').map(decodeURIComponent))
    .reduce((params, [key, value]) => ({ ...params, [key]: value }), {})
}

function params_stringify(params) {
  return Object.keys(params)
    .filter(param => params[param] !== undefined && params[param] !== '')
    .map((param) => [param, params[param]].map(encodeURIComponent).join('='))
    .join('&')
}

function hash_parse(hash_encoded) {
  if (!hash_encoded.startsWith('/')) hash_encoded = `/${hash_encoded}`
  const q = hash_encoded.indexOf('?')
  if (q === -1) {
    return { path: hash_encoded, params: {} }
  } else {
    const [path, params_encoded] = [hash_encoded.slice(0, q), hash_encoded.slice(q+1)]
    const params = params_parse(params_encoded)
    return { path, params }
  }
}

function hash_stringify({ path, params }) {
  const params_encoded = params_stringify(params)
  if (params_encoded === '') {
    return path
  } else {
    return `${path}?${params_encoded}`
  }
}

function hash_get() {
  return (window.location.hash || '#').slice(1)
}

function url_hash_store() {
  let init = hash_get()
  const { path: initPath, params: initParams } = hash_parse(init)
  const { subscribe, update, set } = writable({ path: initPath, params: initParams })
  let lastPath = initPath
  subscribe(({ path, params }) => {
    if (path !== lastPath) {
      // add changes to path to history
      window.location.hash = hash_stringify({ path, params })
      lastPath = path
    } else {
      // don't add changes to params to history, but modify the hash url
      history.replaceState(undefined, undefined, `#${hash_stringify({ path, params })}`)
    }
  })
  window.addEventListener('hashchange', () => set(hash_parse(hash_get())))

  return {
    subscribe,
    update,
    set,
  }
}

const hash = url_hash_store()
export default hash