import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/DGGS.jl/dev/',// TODO: replace this in makedocs!
  title: 'DGGS.jl',
  description: 'Documentation for DGGS.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }]],
  ignoreDeadLinks: true,

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
        md.use(mathjax3),
        md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    }
  },
  themeConfig: {
    outline: 'deep',
    logo: {
      src: 'icon.drawio.svg',
      width: 40,
      height: 40
    },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav: [
{ text: 'Home', link: '/index' },
{ text: 'Guide', collapsed: false, items: [
{ text: 'Get Started', link: '/get_started' },
{ text: 'Background', link: '/background' },
{ text: 'Convert', link: '/convert' },
{ text: 'Select', link: '/select' },
{ text: 'Plot', link: '/plot' }]
 },
{ text: 'API', link: '/api' }
]
,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Guide', collapsed: false, items: [
{ text: 'Get Started', link: '/get_started' },
{ text: 'Background', link: '/background' },
{ text: 'Convert', link: '/convert' },
{ text: 'Select', link: '/select' },
{ text: 'Plot', link: '/plot' }]
 },
{ text: 'API', link: '/api' }
]
,
    editLink: { pattern: "https://github.com/danlooo/DGGS.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/danlooo/DGGS.jl' }
    ],
    footer: {
      message: '<a href="https://www.bgc-jena.mpg.de/en"><img src="logo-mpi-bgc.svg" class = "footer-logo"/></a> <a href="https://earthmonitor.org/"><img src="logo-open-earth-monitor.png" class = "footer-logo"/></a> <a href="https://cordis.europa.eu/project/id/101059548"><img src="logo-eu.png" class = "footer-logo"/></a>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
