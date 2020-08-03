import { writable } from 'svelte/store'

function url_hash_store() {
  const init = (window.location.hash || '#').slice(1)
  const { subscribe, update, set } = writable(init)
  subscribe((value) => window.location.hash = value)
  window.addEventListener('hashchange', () => set((window.location.hash || '#').slice(1)))

  return {
    subscribe,
    update,
    set,
  }
}

export const hash = url_hash_store()
