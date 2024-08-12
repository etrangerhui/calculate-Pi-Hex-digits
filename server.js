// server.js
const express = require('express');
const WebSocket = require('ws');
const path = require('path');
	var perc=0.,per=0.; //����pi�ڼ�λ�Ľ���%
	var fperc=0.,fper;  //������͹�ʽ�Ľ���% , �õ���ʱ�����,���ڼ��ʱ�������
	var epsf=1000;   //ʱ����1��
	var testbig=0000;
	var u=1;
	var RADIX=3;
	var RAD=(1<<RADIX);
	var PS=10;           //ģ��16����С�����7λ
	var NHX=1;          //�ɿ�����ֻ��8λ,!!!��һ����FFFFFFFF������,����Ҳ�ǲ��ɿ���
	var MAX=(u<<((PS-6)<<2))*(u<<((NHX+6)<<2));  //����=16^15
const app = express();
let count = 0; 

let totalSum = 0;
let clients = []; // ���ڱ���ͻ�������


// Serve static files (e.g., client.html)
app.use(express.static(path.join(__dirname, 'public')));

// Create WebSocket server
const wss = new WebSocket.Server({ noServer: true });
var result0='3.',result='';
// Handle WebSocket connection
wss.on('connection', (ws) => {
    const myid = clients.length;
    clients.push(ws);

    ws.send(JSON.stringify({ type: 'register', myid }));

    ws.on('message', (message) => {
        const data = JSON.parse(message);
        if (data.type === 'compute') {
            totalSum += data.pid;
            totalSum %= MAX;
            if (totalSum < 0) totalSum += MAX;

            result= totalSum.toString(16);
            if(result.length<PS+NHX)result="0";
            else result=result[0];
            //ws.send(JSON.stringify({ type: 'result', result }));
            
        }
    });

    ws.on('close', () => {
        clients = clients.filter(client => client !== ws);
    });
});

// Handle HTTP upgrade requests for WebSocket
app.server = app.listen(3000);
app.server.on('upgrade', (request, socket, head) => {
    wss.handleUpgrade(request, socket, head, (ws) => {
        wss.emit('connection', ws, request);
    });
});
app.use((req, res, next) => {  
  res.header('Access-Control-Allow-Origin', '*');  
  res.header('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE');  
  res.header('Access-Control-Allow-Headers', 'Content-Type');  
  next();  
});
// Handle the "calculate" request
app.get('/calculate', (req, res) => {
    testbig += 1;totalSum=0;

    // Broadcast numprocs and testbig to all clients
     clients.forEach(client => {
        client.send(JSON.stringify({ type: 'calculate', numprocs: clients.length, testbig }));
    });

    res.json(result);
    result0+=result;
    count++;  
    
    // ����Ƿ�ﵽ��80���ַ��ı�������ģ80����0��  
    if (count % 80 === 0) {  
        // ����ǣ��������ǰ�ۻ���result��������result��count  
        console.log(result0);  
        result0 = ''; // ����result��Ϊ��һ����׼��  
        count = 0; // ���ü�����  
    }  
});

console.log('Server running on http://localhost:3000');
