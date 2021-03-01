
// https://stackoverflow.com/questions/3426404/create-a-hexadecimal-colour-based-on-a-string-with-javascript
export function hashCode(str) {
  var hash = 0;
  for (var i = 0;i < str.length;i++) {
    hash = str.charCodeAt(i) + ((hash << 5) - hash);
  }
  return hash;
}

export function intToRGB(i) {
  var c = (i & 0x00FFFFFF)
    .toString(16)
    .toUpperCase();
  return '#' + ("00000".substring(0, 6 - c.length) + c);
}

export function localize_appyter_image(base_url, { name, image }) {
  if (/^https?:\/\//.exec(image) !== null) {
    return image
  } else {
    return `${base_url}/${name}/static/${image}`
  }
}

export function set_gte(A, B) {
  for (const b of B) {
    if (!A.has(b)) return false
  }
  return true
}
