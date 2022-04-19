<#macro mainLayout active bodyClass>
<!doctype html>
<html>
<head>
    <meta charset="utf-8">
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <meta name="robots" content="noindex, nofollow">
    <meta name="viewport" content="width=device-width,initial-scale=1">

    <title>${msg("accountManagementTitle")}</title>
    <link rel="icon" href="/favicon.ico">
    <#if properties.stylesCommon?has_content>
        <#list properties.stylesCommon?split(' ') as style>
            <link href="${url.resourcesCommonPath}/${style}" rel="stylesheet" />
        </#list>
    </#if>
    <#if properties.styles?has_content>
        <#list properties.styles?split(' ') as style>
            <link href="${url.resourcesPath}/${style}" rel="stylesheet" />
        </#list>
    </#if>
    <#if properties.scripts?has_content>
        <#list properties.scripts?split(' ') as script>
            <script type="text/javascript" src="${url.resourcesPath}/${script}"></script>
        </#list>
    </#if>

    <style>
      /* sticky footer */
      html {
        margin: 0 !important;
        padding: 0 !important;
      }
      body {
        margin: 0 !important;
        padding: 0 !important;
        display: flex;
        min-height: 100vh;
        min-width: 540px;
        flex-direction: column;
        font-size: 16px;
        font-weight: 400;
        line-height: 1.5;
        background-color: #f5f5f5;
      }
      
      .bg-white {
        background-color: white!important;
      }

      .flex-grow {
        flex: 1 0 auto;
      }

      .flex-grow-1 {
        flex: 1!important;
      }

      .img-fluid {
        max-width: 100%;
        height: auto;
      }

      .p-4 {
        padding: 1.5rem!important;
      }
      .pb-2 {
        padding-bottom: 0.5rem!important;
      }
      .pb-3 {
        padding-bottom: 1.0rem!important;
      }
      .pt-4 {
        padding-top: 1.5rem!important;
      }
      .mb-4 {
        margin-bottom: 1.5rem!important;
      }
      .mt-4 {
        margin-top: 1.5rem!important;
      }
      .mx-0 {
        margin-left: 0!important;
        margin-right: 0!important;
      }
      .my-2 {
        margin-top: 0.5rem!important;
        margin-bottom: 0.5rem!important;
      }
      .my-4 {
        margin-top: 1.5rem!important;
        margin-bottom: 1.5rem!important;
      }
      .justify-content-center {
        justify-content: center!important;
      }
      .align-content-center {
        align-content: center!important;
      }
      .align-items-center {
        align-items: center!important;
      }
      .d-flex {
        display: flex;
      }
      .d-flex-column {
        display: flex;
        flex-direction: column;
      }
      .d-inline-block {
        display: inline-block;
      }
      .text-center: {
        text-align: center!important;
      }
      .text-left: {
        text-align: left!important;
      }
      
      .container-fluid {
        margin-left: 0!important;
        margin-right: 0!important;
      }
    
      /* header */
      span.text-muted > a {
        text-decoration: none;
        color: inherit;
      }
      span.text-muted > a:hover {
        text-decoration: underline;
        cursor: pointer;
      }
    </style>
</head>

<body class="admin-console user ${bodyClass}">
    <div class="container-fluid bg-white pb-2 mb-4">
      <div class="container">
        <div class="row">
          <div class="col-md-12 col-lg-4">
            <h1 class="m-0">
              <a href="/#/">
                <img
                  src="/images/appyters_logo.svg"
                  class="img-responsive w-100 p-4"
                  alt="Appyters"
                />
              </a>
            </h1>
          </div>
          <div class="col-sm-12 col-lg-8 offset-xl-1 col-xl-7 text-center my-4">
            <span class="text-muted">
              <a
                href="/#/what-is-an-appyter/"
              >
                What is an Appyter?
              </a>
                |
              <a
                href="/#/creating-appyters/"
              >
                Creating Appyters
              </a>
              | 
              <a
                href="/#/publishing-appyters/"
              >
                Publishing Appyters
              </a>
              | 
              <a
                href="/#/about/"
              >
                About
              </a>
              | 
              <a
                href="/#/account/"
              >
                Account
              </a>
              <#if realm.internationalizationEnabled>
                |
                <div class="kc-dropdown" id="kc-locale-dropdown">
                  <a href="#" id="kc-current-locale-link">${locale.current}</a>
                  <ul>
                    <#list locale.supported as l>
                      <li class="kc-dropdown-item"><a href="${l.url}">${l.label}</a></li>
                    </#list>
                  </ul>
                </div>
              </#if>
              <#if referrer?has_content && referrer.url?has_content>
                |
                <a href="${referrer.url}" id="referrer" style="font-weight: 600">${msg("backTo",referrer.name)}</a>
              </#if>
            </span>
          </div>
        </div>
      </div>
    </div>
    <div class="container-fluid flex-grow">
      
      <div class="container">
          <div class="bs-sidebar col-sm-3">
              <ul>
                  <li class="<#if active=='account'>active</#if>"><a href="${url.accountUrl}">${msg("account")}</a></li>
                  <#if features.passwordUpdateSupported><li class="<#if active=='password'>active</#if>"><a href="${url.passwordUrl}">${msg("password")}</a></li></#if>
                  <li class="<#if active=='totp'>active</#if>"><a href="${url.totpUrl}">${msg("authenticator")}</a></li>
                  <#if features.identityFederation><li class="<#if active=='social'>active</#if>"><a href="${url.socialUrl}">${msg("federatedIdentity")}</a></li></#if>
                  <li class="<#if active=='sessions'>active</#if>"><a href="${url.sessionsUrl}">${msg("sessions")}</a></li>
                  <li class="<#if active=='applications'>active</#if>"><a href="${url.applicationsUrl}">${msg("applications")}</a></li>
                  <#if features.log><li class="<#if active=='log'>active</#if>"><a href="${url.logUrl}">${msg("log")}</a></li></#if>
                  <#if realm.userManagedAccessAllowed && features.authorization><li class="<#if active=='authorization'>active</#if>"><a href="${url.resourceUrl}">${msg("myResources")}</a></li></#if>
              </ul>
          </div>

          <div class="col-sm-9 content-area">
              <#if message?has_content>
                  <div class="alert alert-${message.type}">
                      <#if message.type=='success' ><span class="pficon pficon-ok"></span></#if>
                      <#if message.type=='error' ><span class="pficon pficon-error-circle-o"></span></#if>
                      <span class="kc-feedback-text">${kcSanitize(message.summary)?no_esc}</span>
                  </div>
              </#if>

              <#nested "content">
          </div>
      </div>
  </div>
  <div class="footer mt-4 pt-4 bg-white">
    <div class="d-flex col-md-5 col-sm-12 justify-content-center align-items-center">
      <p class="d-inline-block text-left">
        <a style="color: #555;" href="mailto:avi.maayan@mssm.edu">Contact Us</a><br />
        <a style="color: #555;" href="https://github.com/MaayanLab/appyter-catalog/blob/main/LICENSE">Usage License</a><br />
        <a style="color: #555;" href="/#/what-is-an-appyter/">Appyter Documentation</a><br />
      </p>
    </div>
    <div class="col-md-2 col-xs-3 text-center">
      <a href="https://icahn.mssm.edu/research/bioinformatics" target="_blank">
        <img class="rounded" src="/images/icahn_cb.png" style="height: 5rem;">
      </a>
    </div>
    <div class="col-md-2 col-xs-3 text-center">
      <a href="https://labs.icahn.mssm.edu/maayanlab/" target="_blank">
        <img class="rounded" src="/images/maayanlab_logo.png" style="height: 5rem;">
      </a>
    </div>
    <div class="d-flex-column col-md-2 col-xs-3 text-center">
      <div class="my-2">
        <a class="badge badge-secondary px-2" href="https://github.com/MaayanLab/appyter-catalog" target="_blank">
          <span class="badge badge-light">
            <img src="/images/GitHub-Mark.png" style="width: 1rem;">
          </span>
          View source code
        </a>
      </div>
      <div class="my-2">
        <a class="badge badge-secondary px-2" href="https://github.com/MaayanLab/appyter-catalog/issues/new" target="_blank">
          <span class="badge badge-light">
            <img src="/images/GitHub-Mark.png" style="width: 1rem;">
          </span>
          Submit an issue
        </a>
      </div>
    </div>
  </div>

</body>
</html>
</#macro>