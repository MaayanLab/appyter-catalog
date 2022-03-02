<script>
  import clipboardIcon from 'bootstrap-icons/icons/clipboard.svg'

  export let code = ''

  let copy_feedback = false

  let code_lines
  $: {
    let left_spaces
    let lines = code.split('\n')
    let fixed_lines = []
    for (const line of lines) {
      if (left_spaces === undefined) {
        const line_trimmed = line.replace(/^ */g, '')
        // remove blank lines at the beginning
        if (line_trimmed === '') continue
        // figure out indent for dedenting
        left_spaces = line.length - line_trimmed.length
      }
      const fixed_line = line.slice(left_spaces).replace(/ *$/g, '')
      fixed_lines.push(fixed_line)
    }
    // remove blank lines at the end
    while (fixed_lines.length > 0 && fixed_lines[fixed_lines.length - 1] === '') {
      fixed_lines.pop()
    }
    code_lines = fixed_lines
  }
</script>

<div
  class="row mx-0 my-3"
  style="background-color: #fafafa; "
>
  <div class="col p-3" style="overflow-x: auto;">
    <code style="white-space: nowrap">
      {#each code_lines as line}
        {#if line.startsWith('#')}
          <span style="color: green">{line}</span>
        {:else}
          <span>{line}</span>
        {/if}
        <br />
      {/each}
    </code>
  </div>
  <div class="p-3">
    <button
      class="btn btn-sm"
      class:btn-success={copy_feedback}
      title="Copy to clipboard"
      on:click={() => {
        navigator.clipboard.writeText(code_lines.join('\n'))
          .then(() => {
            copy_feedback = true
            setTimeout(() => copy_feedback = false, 1000)
          })
      }}
    ><img src={clipboardIcon} alt="Clipboard" /></button>
  </div>
</div>