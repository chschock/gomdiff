# Introduction

This is a golang port of the sparse linear assignment problem solver in my project git@github.com:chschock/mdiff.git .

The interface is a webservice api that. The frontend of mdiff can also send requests to gomdiff and therefore provides a usage example.

You have to install a browser plugin to allow cross-origin resource sharing, if you want to run in on localhost without webserver.

# Performance

The performance is significantly faster (3x) than javascript, although I expected more.

The time statistics logged by the server do not include communication time. The packing and unpacking of the json response is actually quite slow, but for the current version of the frontend this data is necessary.
