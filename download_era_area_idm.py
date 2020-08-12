# -*- coding: utf-8 -*-
#"area": [30.5,113.5,24,118.5]
#year
#dd

from __future__ import absolute_import, division, print_function, unicode_literals
from subprocess import call
import sys
import json
import time
import os
import logging
import uuid
import requests

try:
    from urllib.parse import urljoin
except ImportError:
    from urlparse import urljoin



def bytes_to_string(n):
    u = ['', 'K', 'M', 'G', 'T', 'P']
    i = 0
    while n >= 1024:
        n /= 1024.0
        i += 1
    return '%g%s' % (int(n * 10 + 0.5) / 10.0, u[i])


def read_config(path):
    config = {}
    with open(path) as f:
        for l in f.readlines():
            if ':' in l:
                k, v = l.strip().split(':', 1)
                if k in ('url', 'key', 'verify'):
                    config[k] = v.strip()
    return config


def toJSON(obj):

    to_json = getattr(obj, "toJSON", None)
    if callable(to_json):
        return to_json()

    if isinstance(obj, (list, tuple)):
        return [toJSON(x) for x in obj]

    if isinstance(obj, dict):
        r = {}
        for k, v in obj.items():
            r[k] = toJSON(v)
        return r

    return obj


class Result(object):

    def __init__(self, client, reply):

        self.reply = reply

        self._url = client.url

        self.session = client.session
        self.robust = client.robust
        self.verify = client.verify
        self.cleanup = client.delete

        self.debug = client.debug
        self.info = client.info
        self.warning = client.warning
        self.error = client.error
        self.sleep_max = client.sleep_max
        self.retry_max = client.retry_max

        self.timeout = client.timeout
        self.progress = client.progress

        self._deleted = False

    def toJSON(self):
        r = dict(resultType='url',
                 contentType=self.content_type,
                 contentLength=self.content_length,
                 location=self.location)
        return r

    def _download(self, url, size, target):

        if target is None:
            target = url.split('/')[-1]

        self.info("Downloading %s to %s (%s)", url, target, bytes_to_string(size))
        fpath=os.path.split(target)        
        #call(['IDM.exe', '/d',url, '/p',fpath[0], '/f', fpath[1], '/n', '/s'])
        os.system('IDM.exe /d '+url+' /p '+fpath[0]+' /f '+fpath[1]+' /n /s')
        
        
        return target

    def download(self, target=None):
        return self._download(self.location,
                              self.content_length,
                              target)

    @property
    def content_length(self):
        return int(self.reply['content_length'])

    @property
    def location(self):
        return urljoin(self._url, self.reply['location'])

    @property
    def content_type(self):
        return self.reply['content_type']

    def __repr__(self):
        return "Result(content_length=%s,content_type=%s,location=%s)" % (self.content_length,
                                                                          self.content_type,
                                                                          self.location)

    def check(self):
        self.debug("HEAD %s", self.location)
        metadata = self.robust(self.session.head)(self.location,
                                                  verify=self.verify,
                                                  timeout=self.timeout)
        metadata.raise_for_status()
        self.debug(metadata.headers)
        return metadata

    def delete(self):

        if self._deleted:
            return

        if 'request_id' in self.reply:
            rid = self.reply['request_id']

            task_url = '%s/tasks/%s' % (self._url, rid)
            self.debug("DELETE %s", task_url)

            delete = self.session.delete(task_url, verify=self.verify)
            self.debug("DELETE returns %s %s", delete.status_code, delete.reason)

            try:
                delete.raise_for_status()
            except Exception:
                self.warning("DELETE %s returns %s %s",
                             task_url, delete.status_code, delete.reason)

            self._deleted = True

    def __del__(self):
        try:
            if self.cleanup:
                self.delete()
        except Exception as e:
            print(e)


class Client(object):

    logger = logging.getLogger('cdsapi')

    def __init__(self,
                 url=os.environ.get('CDSAPI_URL'),
                 key=os.environ.get('CDSAPI_KEY'),
                 quiet=False,
                 debug=False,
                 verify=None,
                 timeout=60,
                 progress=True,
                 full_stack=False,
                 delete=True,
                 retry_max=500,
                 sleep_max=120,
                 info_callback=None,
                 warning_callback=None,
                 error_callback=None,
                 debug_callback=None,
                 ):

        if not quiet:

            if debug:
                level = logging.DEBUG
            else:
                level = logging.INFO

            logging.basicConfig(level=level,
                                format='%(asctime)s %(levelname)s %(message)s')

        dotrc = os.environ.get('CDSAPI_RC', os.path.expanduser('~/.cdsapirc'))

        if url is None or key is None:
            if os.path.exists(dotrc):
                config = read_config(dotrc)

                if key is None:
                    key = config.get('key')

                if url is None:
                    url = config.get('url')

                if verify is None:
                    verify = int(config.get('verify', 1))

        if url is None or key is None or key is None:
            raise Exception('Missing/incomplete configuration file: %s' % (dotrc))

        self.url = url
        self.key = key

        self.quiet = quiet
        self.progress = progress and not quiet

        self.verify = True if verify else False
        self.timeout = timeout
        self.sleep_max = sleep_max
        self.retry_max = retry_max
        self.full_stack = full_stack
        self.delete = delete
        self.last_state = None

        self.debug_callback = debug_callback
        self.warning_callback = warning_callback
        self.info_callback = info_callback
        self.error_callback = error_callback

        self.session = requests.Session()
        self.session.auth = tuple(self.key.split(':', 2))

        self.debug("CDSAPI %s", dict(url=self.url,
                                     key=self.key,
                                     quiet=self.quiet,
                                     verify=self.verify,
                                     timeout=self.timeout,
                                     progress=self.progress,
                                     sleep_max=self.sleep_max,
                                     retry_max=self.retry_max,
                                     full_stack=self.full_stack,
                                     delete=self.delete
                                     ))

    def retrieve(self, name, request, target=None):
        result = self._api('%s/resources/%s' % (self.url, name), request, 'POST')
        if target is not None:
            result.download(target)
        return result

    def service(self, name, *args, **kwargs):
        self.delete = False  # Don't delete results
        name = '/'.join(name.split('.'))
        request = toJSON(dict(args=args, kwargs=kwargs))
        result = self._api('%s/tasks/services/%s/clientid-%s' % (self.url, name, uuid.uuid4().hex), request, 'PUT')
        return result

    def workflow(self, code, *args, **kwargs):
        params = dict(code=code,
                      args=args,
                      kwargs=kwargs,
                      workflow_name='application')
        return self.service("tool.toolbox.orchestrator.run_workflow", params)

    def identity(self):
        return self._api('%s/resources' % (self.url,), {})

    def _api(self, url, request, method):

        session = self.session

        self.info("Sending request to %s", url)
        self.debug("%s %s %s", method, url, json.dumps(request))

        if method == 'PUT':
            action = session.put
        else:
            action = session.post

        result = self.robust(action)(url,
                                     json=request,
                                     verify=self.verify,
                                     timeout=self.timeout)
        reply = None

        try:
            result.raise_for_status()
            reply = result.json()
        except Exception:

            if reply is None:
                try:
                    reply = result.json()
                except Exception:
                    reply = dict(message=result.text)

            self.debug(json.dumps(reply))

            if 'message' in reply:
                error = reply['message']

                if 'context' in reply and 'required_terms' in reply['context']:
                    e = [error]
                    for t in reply['context']['required_terms']:
                        e.append("To access this resource, you first need to accept the terms"
                                 "of '%s' at %s" % (t['title'], t['url']))
                    error = '. '.join(e)
                raise Exception(error)
            else:
                raise

        sleep = 1

        while True:

            self.debug("REPLY %s", reply)

            if reply['state'] != self.last_state:
                self.info("Request is %s" % (reply['state'],))
                self.last_state = reply['state']

            if reply['state'] == 'completed':
                self.debug("Done")

                if 'result' in reply:
                    return reply['result']

                return Result(self, reply)

            if reply['state'] in ('queued', 'running'):
                rid = reply['request_id']

                self.debug("Request ID is %s, sleep %s", rid, sleep)
                time.sleep(sleep)
                sleep *= 1.5
                if sleep > self.sleep_max:
                    sleep = self.sleep_max

                task_url = '%s/tasks/%s' % (self.url, rid)
                self.debug("GET %s", task_url)

                result = self.robust(session.get)(task_url,
                                                  verify=self.verify,
                                                  timeout=self.timeout)
                result.raise_for_status()
                reply = result.json()
                continue

            if reply['state'] in ('failed',):
                self.error("Message: %s", reply['error'].get('message'))
                self.error("Reason:  %s", reply['error'].get('reason'))
                for n in reply.get('error', {}).get('context', {}).get('traceback', '').split('\n'):
                    if n.strip() == '' and not self.full_stack:
                        break
                    self.error("  %s", n)
                raise Exception("%s. %s." % (reply['error'].get('message'), reply['error'].get('reason')))

            raise Exception('Unknown API state [%s]' % (reply['state'],))

    def info(self, *args, **kwargs):
        if self.info_callback:
            self.info_callback(*args, **kwargs)
        else:
            self.logger.info(*args, **kwargs)

    def warning(self, *args, **kwargs):
        if self.warning_callback:
            self.warning_callback(*args, **kwargs)
        else:
            self.logger.warning(*args, **kwargs)

    def error(self, *args, **kwargs):
        if self.error_callback:
            self.error_callback(*args, **kwargs)
        else:
            self.logger.error(*args, **kwargs)

    def debug(self, *args, **kwargs):
        if self.debug_callback:
            self.debug_callback(*args, **kwargs)
        else:
            self.logger.debug(*args, **kwargs)

    def _download(self, results, targets=None):

        if isinstance(results, Result):
            if targets:
                path = targets.pop(0)
            else:
                path = None
            return results.download(path)

        if isinstance(results, (list, tuple)):
            return [self._download(x, targets) for x in results]

        if isinstance(results, dict):

            if 'location' in results and 'contentLength' in results:
                reply = dict(location=results['location'],
                             content_length=results['contentLength'],
                             content_type=results.get('contentType'))

                if targets:
                    path = targets.pop(0)
                else:
                    path = None

                return Result(self, reply).download(path)

            r = {}
            for k, v in results.items():
                r[v] = self._download(v, targets)
            return r

        return results

    def download(self, results, targets=None):
        if targets:
            # Make a copy
            targets = [t for t in targets]
        return self._download(results, targets)

    def remote(self, url):
        r = requests.head(url)
        reply = dict(location=url,
                     content_length=r.headers['Content-Length'],
                     content_type=r.headers['Content-Type'])
        return Result(self, reply)

    def robust(self, call):

        def retriable(code, reason):

            if code in [requests.codes.internal_server_error,
                        requests.codes.bad_gateway,
                        requests.codes.service_unavailable,
                        requests.codes.gateway_timeout,
                        requests.codes.too_many_requests,
                        requests.codes.request_timeout]:
                return True

            return False

        def wrapped(*args, **kwargs):
            tries = 0
            while tries < self.retry_max:
                try:
                    r = call(*args, **kwargs)
                except requests.exceptions.ConnectionError as e:
                    r = None
                    self.warning("Recovering from connection error [%s], attemps %s of %s",
                                 e, tries, self.retry_max)

                if r is not None:
                    if not retriable(r.status_code, r.reason):
                        return r
                    self.warning("Recovering from HTTP error [%s %s], attemps %s of %s",
                                 r.status_code, r.reason, tries, self.retry_max)

                tries += 1

                self.warning("Retrying in %s seconds", self.sleep_max)
                time.sleep(self.sleep_max)
                self.info("Retrying now...")

        return wrapped


if __name__ == '__main__':     
    argc = len(sys.argv)  
    if argc<9:
        print ('Need to input year mbegin mend latmin lonmin latmax lonmax outdir israin')
        sys.exit(1) 
        
    year=sys.argv[1] 
    mbegin=int(sys.argv[2])
    mend=int(sys.argv[3])
    latmin=float(sys.argv[4])
    lonmin=float(sys.argv[5])
    latmax=float(sys.argv[6])
    lonmax=float(sys.argv[7])
    outdir=sys.argv[8]
    israin=int(sys.argv[9])


    c = Client()  
    for dd in range(mbegin,mend+1):
        
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type':'reanalysis',
                'format':'netcdf',
                'variable':[
                    'specific_cloud_ice_water_content','specific_cloud_liquid_water_content',
                    'specific_rain_water_content','specific_snow_water_content'
                ],
                'pressure_level':[
                    '100','150',
                    '200','250','300',
                    '350','400','450',
                    '500','550','600',
                    '650','700','750',
                    '775','800','825',
                    '850','875','900',
                    '925','950','975',
                    '1000'
                ],
                'year':year,
                'month':'%02d'%(dd),
                'day':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12',
                '13','14','15',
                '16','17','18',
                '19','20','21',
                '22','23','24',
                '25','26','27',
                '28','29','30',
                '31'
                ],
                'time':[
                    '00:00','01:00','02:00',
                    '03:00','04:00','05:00',
                    '06:00','07:00','08:00',
                    '09:00','10:00','11:00',
                    '12:00','13:00','14:00',
                    '15:00','16:00','17:00',
                    '18:00','19:00','20:00',
                    '21:00','22:00','23:00'
                ],
                "area": [latmax,lonmin,latmin,lonmax],
                
            },
            outdir+'\\era5_'+year+'%02d'%(dd)+'_qc.nc')
        
        #######
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type':'reanalysis',
                'format':'netcdf',
                'variable':[
                    'specific_humidity',
                    'u_component_of_wind','v_component_of_wind'
                ],
                'pressure_level':[
                    '100','150',
                    '200','250','300',
                    '350','400','450',
                    '500','550','600',
                    '650','700','750',
                    '775','800','825',
                    '850','875','900',
                    '925','950','975',
                    '1000'
                ],
                'year':year,
                'month':'%02d'%(dd),
                'day':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12',
                '13','14','15',
                '16','17','18',
                '19','20','21',
                '22','23','24',
                '25','26','27',
                '28','29','30',
                '31'
                ],
                'time':[
                    '00:00','01:00','02:00',
                    '03:00','04:00','05:00',
                    '06:00','07:00','08:00',
                    '09:00','10:00','11:00',
                    '12:00','13:00','14:00',
                    '15:00','16:00','17:00',
                    '18:00','19:00','20:00',
                    '21:00','22:00','23:00'
                ],
                "area": [latmax,lonmin,latmin,lonmax],
                
            },
            outdir+'\\era5_'+year+'%02d'%(dd)+'_quv.nc')
           
        #######
        if israin==1:
            c.retrieve(
                    'reanalysis-era5-land-monthly-means',
                    {
                        'format':'netcdf',
                        'variable':[
                            'total_evaporation','snow_evaporation','total_precipitation'
                        ],
                        'year':year,
                        'month':'%02d'%(dd),
                        "area": [latmax,lonmin,latmin,lonmax],
                        'time':'00:00',
                        'product_type':'monthly_averaged_reanalysis'
                    },
                    outdir+'\\era5_'+year+'%02d'%(dd)+'_rain.nc')
        