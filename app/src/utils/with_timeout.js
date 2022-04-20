export default function with_timeout(promise, timeout) {
  return (...args) => new Promise((resolve, reject) => {
    let ctx = { done: false }
    setTimeout(() => {
      if (ctx.done) return
      ctx.done = true
      reject(new Error('Timeout'))
    }, timeout)
    promise(...args)
      .then((value) => {
        if (ctx.done) return
        ctx.done = true
        resolve(value)
      })
      .catch((e) => {
        if (ctx.done) return
        ctx.done = true
        reject(e)
      })
  })
}