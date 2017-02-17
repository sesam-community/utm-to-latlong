===========================
UTM to LatLong microservice
===========================

A python micro service template for transforming a JSON entity stream. This service is designed to be used as a microservice system with
the `HTTP transform <https://docs.sesam.io/configuration.html#the-http-transform>`_ in a Sesam service instance.

.. image:: https://travis-ci.org/sesam-community/utmtolatlong-microservice.svg?branch=master
   :alt: Build Status
   :target: https://travis-ci.org/sesam-community/utmtolatlong-microservice

It converts UTM33 (EUREF89 aka GAB) coordinates to LatLong (GPS coordinates aka WGS84).


Running locally in a virtual environment
----------------------------------------

::

  cd utmtolatlong-microservice/service
  virtualenv --python=python3 venv
  . venv/bin/activate
  pip install -r requirements.txt

  python transform-service.py
   * Running on http://0.0.0.0:5001/ (Press CTRL+C to quit)
   * Restarting with stat
   * Debugger is active!
   * Debugger pin code: 260-787-156

The service listens on port 5001 on localhost.

Running in Docker
-----------------

::

  cd utmtolatlong-microservice
  docker build -t utmtolatlong-microservice .
  docker run -it --rm --name utmtolatlong-microservice -p 5001:5001 utmtolatlong-microservice

Get the IP from docker:

::

  docker inspect -f '{{.Name}} - {{.NetworkSettings.IPAddress }}' utmtolatlong-microservice

Example
-------
  
JSON entities can be posted to 'http://localhost:5001/transform'. The result is streamed back to the client. Exchange "localhost" with the Docker IP if running in Docker.

::

   $ curl -s -XPOST 'http://localhost:5001/transform' -H "Content-type: application/json" -d '[{ "_id": "jane", "northing": "12344", "easting": "6543", "zone": "32"}]' | jq -S .
   [
     {
       "_id": "jane",
       "message": "Hello world!",
       "name": "Jane Doe",
       "northing": "12344",
       "easting": "6543",
       "zone": "32",
       "lat": 0.11134243423572525,
       "lon": 4.569868852874141
     }
   ]

Note the example uses `curl <https://curl.haxx.se/>`_ to send the request and `jq <https://stedolan.github.io/jq/>`_ prettify the response.

Configuration
-------------

You can configure the service with the following environment variables:

=====================  =====================================================================================   ==========
Variable               Description                                                                             Default


``EASTING_PROPERTY``   The name of the property holding the "easting" value. The value itself can be
                       either a string, float, int or decimal. If a string, it must be castable to a float.    "easting"

``NORTHING_PROPERTY``  The name of the property holding the "northing" value. The value itself can be either   "northing"
                       a string, float, int or decimal. If a string, it must be castable to a float.          

``ZONE_PROPERTY``      The name of the property holding the "zone" value. The value itself can be either a     "zone"
                       string, float, int or decimal. It must be castable to an int.

``ZONE_DEFAULT``       The default "zone" value if none are found in the entity. The value itself can be       "32"
                       either a string, float, int or decimal. It must be castable to an int. 

``HEMI_PROPERTY``      The name of the property holding the "hemisphere" value. The value itself can be        "hemi"
                       either a string, float, int or decimal. It must be castable to an int. A 0 value
                       means northern hemisphere, any other value means southern hemisphere.

``HEMI_DEFAULT``       The default "hemisphere" value. The value itself can be either a string, float,         "0"
                       int or decimal. It must be castable to an int. A 0 value means northern hemisphere,
                       any other value means southern hemisphere.

``LAT_PROPERTY``       The name of the property that will hold the converted "latitude" value.                 "lat" 

``LONG_PROPERTY``      The name of the property that will hold the converted "longitude" value.                "long"

``INCLUDE_LAT_LONG``   A flag to indicate wheter to include a lat_long value in the returned entity.           "false"
                       A value of "true" means yes. Any other value means no.

``LAT_LONG_PROPERTY``  The name of the property that will hold the converted "lat_long" value                  "lat_long"
                       ("lat, long"). It is only inserted in the returned entity if ``INCLUDE_LAT_LONG``
                       evaluates to yes.
=====================  =====================================================================================   ==========

When running in Docker you can either specify this in a file (see https://docs.docker.com/compose/env-file/) or on the command line with "docker run .. -e VAR1=VAL1 -e VAR2=VAL2 .."
