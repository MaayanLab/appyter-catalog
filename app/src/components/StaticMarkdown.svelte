<script>
  export let data
</script>

<style>
  :global(.md > ol > li) {
    list-style-type: decimal;
  }
  :global(.md > ol > li > ol > li) {
    list-style-type: lower-roman;
  }
  :global(.md > ol > li > ol > li > ol > li) {
    list-style-type: lower-alpha;
  }
  :global(.md p) {
    text-align: justify;
  }
</style>

<div class="md">
  {@html data}
</div>
