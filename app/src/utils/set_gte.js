export default function set_gte(A, B) {
  for (const b of B) {
    if (!A.has(b)) return false
  }
  return true
}
