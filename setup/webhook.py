import json
import requests

f = open("tokens/NGROK.TOKEN",'r')
ngrok_token = f.read().rstrip('\n')
f.close()

f = open("tokens/GIT.TOKEN",'r')
git_token = f.read().rstrip('\n')
f.close()

f = open("tokens/JENKINS.TOKEN",'r')
jenkins_token = f.read().rstrip('\n')
f.close()

headers={"Authorization": "Bearer " + str(ngrok_token), "Ngrok-version": "2"}
endpoints = requests.get('https://api.ngrok.com/endpoints', headers=headers)
url_ngrok=endpoints.json()['endpoints'][0]['public_url']
url_ngrok=url_ngrok.split('/')[2]

url="https://matthew:"+jenkins_token+"@"+url_ngrok+"/github-webhook/"
hook_id="405325584"
git_headers={"Accept": "application/vnd.github+json", "Authorization": "Bearer "+git_token, "X-GitHub-Api-Version": "2022-11-28"}
#git_params={"name":"web","active":True,"events":["push","pull_request"],"config":{"url":url,"content_type":"json","insecure_ssl":"0"}}
#requests.post("https://api.github.com/repos/09nwalkerm/ROMEG/hooks",headers=git_headers, json=git_params)
git_params={"url":url,"content_type":"json"}
requests.patch("https://api.github.com/repos/09nwalkerm/ROMEG/hooks/"+hook_id+"/config",headers=git_headers, json=git_params)
