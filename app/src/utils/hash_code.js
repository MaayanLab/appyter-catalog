// https://stackoverflow.com/questions/3426404/create-a-hexadecimal-colour-based-on-a-string-with-javascript
export default function hashCode(str) {
  var hash = 0;
  for (var i = 0;i < str.length;i++) {
    hash = str.charCodeAt(i) + ((hash << 5) - hash);
  }
  return hash;
}