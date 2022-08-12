import path from 'path'
import MiniCssExtractPlugin from 'mini-css-extract-plugin'
import HtmlWebpackPlugin from 'html-webpack-plugin'
import webpack from 'webpack'
import TerserWebpackPlugin from 'terser-webpack-plugin'
import indexHtml from './src/index.html.js'

export default function (_env, argv) {
  const __dirname = path.resolve();
  const isProduction = argv.mode === 'production';
  const isDevelopment = !isProduction;

  return {
    devtool: isDevelopment && 'cheap-module-source-map',
    entry: './src/index.js',
    output: {
      path: path.resolve(__dirname, 'dist'),
      filename: 'assets/js/[name].[chunkhash:8].js',
      chunkFilename: 'assets/js/[name].[chunkhash:8].chunk.js',
      assetModuleFilename: 'assets/media/[name].[hash:8][ext]',
      publicPath: '/',
      library: {
        name: 'appyters',
        type: 'umd',
        export: 'default',
      },
    },
    resolve: {
      extensions: ['.mjs', '.js', '.md', '.svelte'],
      mainFields: ['svelte', 'browser', 'module', 'main'],
      alias: {
        '@/public': path.resolve(__dirname, 'public'),
        '@': path.resolve(__dirname, 'src'),
      },
    },
    module: {
      rules: [
        {
          test: /\.svelte$/,
          use: {
            loader: 'svelte-loader',
            options: {
              compilerOptions: {
                dev: isDevelopment,
              },
              emitCss: isProduction,
              hotReload: isDevelopment,
            },
          },
        },
        {
          test: /\.(png|jpg|gif|svg)$/,
          type: 'asset/resource',
        },
        {
          test: /\.md$/,
          use: [
            'html-loader',
            'markdown-loader',
          ],
        },
        {
          test: /\.css$/i,
          use: ["style-loader", "css-loader"],
        },
        {
          // required to prevent errors from Svelte on Webpack 5+
          test: /node_modules\/svelte\/.*\.mjs$/,
          resolve: {
            fullySpecified: false
          }
        },
      ]
    },
    plugins: [
      isProduction && new MiniCssExtractPlugin({
        filename: 'assets/css/[name].css',
        chunkFilename: 'assets/css/[name].[contenthash:8].chunk.css'
      }),
      new HtmlWebpackPlugin({
        templateContent: indexHtml,
        inject: true,
      }),
      new webpack.DefinePlugin({
        'process.env.NODE_ENV': JSON.stringify(
          isProduction ? 'production' : 'development'
        )
      })
    ].filter(Boolean),
    optimization: {
      minimize: isProduction,
      minimizer: [
        new TerserWebpackPlugin({
          terserOptions: {
            compress: {
              comparisons: false
            },
            mangle: {
              safari10: true
            },
            output: {
              comments: false,
              ascii_only: true
            },
            warnings: false
          }
        }),
      ],
      splitChunks: {
        chunks: 'all',
        minSize: 0,
        maxInitialRequests: 10,
        maxAsyncRequests: 10,
        cacheGroups: {
          vendors: {
            test: /[\\/]node_modules[\\/]/,
            name(module, chunks, cacheGroupKey) {
              const packageName = module.context.match(
                /[\\/]node_modules[\\/](.*?)([\\/]|$)/
              )[1];
              return `${cacheGroupKey}.${packageName.replace('@', '')}`;
            }
          },
          common: {
            minChunks: 2,
            priority: -10
          }
        }
      },
      runtimeChunk: 'single'
    },
    devServer: {
      static: {
        directory: path.join(__dirname, 'public'),
      },
      port: 1234,
      compress: true,
      allowedHosts: 'all',
      hot: true,
      client: {
        webSocketURL: 'wss://appyters.u8sand.net/ws',
      },
    }
  };
};