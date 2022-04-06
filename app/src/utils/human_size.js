export default function human_size(size) {
  const units = ['B', 'KB', 'MB', 'GB', 'TB']
  let unit
  for (let i = 0; i < units.length-1; i++) {
    unit = units[i]
    if (size < 1000) {
      break
    } else {
      size /= 1000
    }
  }
  return `${size.toFixed(2)} ${unit}`
}