
export default function intToRGB(i) {
  var c = (i & 0x00FFFFFF)
    .toString(16)
    .toUpperCase();
  return '#' + ("00000".substring(0, 6 - c.length) + c);
}
