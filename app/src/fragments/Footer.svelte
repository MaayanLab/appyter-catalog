<script>
  import icahnLogo from '@/public/images/icahn_cb.png'
  import maayanlabLogo from '@/public/images/maayanlab_logo.png'
  import githubLogo from '@/public/images/GitHub-Mark.png'
</script>

<div class="footer mt-4 pt-4 bg-white">
  <div class="row justify-content-center align-content-center mx-0 pb-3">
    <div class="d-flex col-md-3 col-sm-12 justify-content-center align-items-center">
      <p class="d-inline-block text-left">
        <a style="color: #555; text-decoration: none" href="mailto:avi.maayan@mssm.edu">Contact Us</a><br />
        <a style="color: #555; text-decoration: none" href="https://github.com/MaayanLab/appyter-catalog/blob/main/LICENSE">Usage License</a><br />
        <a style="color: #555; text-decoration: none" href="/#/what-is-an-appyter/">Appyter Documentation</a><br />
      </p>
    </div>
    <div class="col-md-2 col-xs-3 text-center">
      <a href="https://icahn.mssm.edu/research/bioinformatics" target="_blank">
        <img class="rounded" src={icahnLogo} alt="Mount Sinai Logo" style="height: 5rem;">
      </a>
    </div>
    <div class="col-md-2 col-xs-3 text-center">
      <a href="https://labs.icahn.mssm.edu/maayanlab/" target="_blank">
        <img class="rounded" src={maayanlabLogo} alt="Ma'ayanlab Logo" style="height: 5rem;">
      </a>
    </div>
    <div class="d-flex-column col-md-2 col-xs-3 text-center">
      <div class="my-2">
        <a
          class="badge bg-secondary px-2"
          style="text-decoration: none"
          href="https://github.com/MaayanLab/appyter-catalog"
          target="_blank"
        >
          <span class="badge bg-light">
            <img src={githubLogo} alt="Github Logo" style="width: 1rem;">
          </span>
          View source code
        </a>
      </div>
      <div class="my-2">
        <a
          class="badge bg-secondary px-2"
          style="text-decoration: none"
          href="https://github.com/MaayanLab/appyter-catalog/issues/new"
          target="_blank"
        >
          <span class="badge bg-light">
            <img src={githubLogo} alt="Github Logo" style="width: 1rem;">
          </span>
          Submit an issue
        </a>
      </div>
    </div>
  </div>
</div>
