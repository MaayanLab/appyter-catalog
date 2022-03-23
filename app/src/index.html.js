import fs from 'fs'

const appyterJson = JSON.parse(fs.readFileSync('public/appyters.json').toString())

const jsonld = {
  "@context": "http://schema.org",
  "@type":"WebSite",
  "url": "https://appyters.maayanlab.cloud",
  "name": "Appyters",
  "alternativeHeadline": "A catalog of appyter notebooks",
  "description": "Appyters are standalone applications for generating jupyter notebooks, the appyter catalog features many bioinformatics related appyters.",
  "keywords": "Jupyter, Notebook, Machine, Learning, Enrichment, Visualization, Bioinformatics, Maayanlab, Avi Ma'ayan",
  "thumbnailUrl": "https://appyters.maayanlab.cloud/favicon.ico",
  "license": "https://github.com/MaayanLab/appyter-catalog/blob/main/LICENSE",
  "isAccessibleForFree": true,
  "author": {
    "@type": "Organization",
    "name": "Ma'ayan Lab",
    "url": "https://labs.icahn.mssm.edu/maayanlab/",
    "contactPoint": {
      "@type": "ContactPoint",
      "url": "https://labs.icahn.mssm.edu/maayanlab/contact/"
    }
  }
}

export default () => `
<!DOCTYPE html>
<html>
  <head>
    <title>${jsonld.name}</title>
    <meta name="subtitle" content=${JSON.stringify(jsonld.alternativeHeadline)}>
    <meta name="description" content=${JSON.stringify(jsonld.description)}>
    <meta name="keywords" content=${JSON.stringify(jsonld.keywords)}>
    <meta name="author" content=${JSON.stringify(jsonld.author.name)}>
    <meta name="image" content=${JSON.stringify(jsonld.thumbnailUrl)}>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <link rel="icon" href="favicon.ico" sizes="16x16 32x32" type="image/png">
    <script type="application/ld+json">${JSON.stringify(jsonld)}</script>
    ${appyterJson.ga_id ? `
    <script async src="https://www.googletagmanager.com/gtag/js?id=${appyterJson.ga_id}"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', ${JSON.stringify(appyterJson.ga_id)});
    </script>
    ` : ''}
    ${appyterJson.keycloak !== undefined ? `<script src=${JSON.stringify(appyterJson.keycloak.url)}></script>` : ''}
  </head>
  <body>
  </body>
</html>`
