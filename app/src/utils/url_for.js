function params_stringify(params) {
  if (params === undefined) return ''
  return Object.keys(params)
    .filter(param => params[param] !== undefined && params[param] !== null && params[param] !== '')
    .map((param) => [param, params[param]].map(encodeURIComponent).join('='))
    .join('&')
}

export default function url_for({ path, params }) {
  const params_encoded = params_stringify(params)
  if (params_encoded === '') {
    return path
  } else {
    return `${path}?${params_encoded}`
  }
}