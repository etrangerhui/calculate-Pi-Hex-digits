// server.js
const express = require('express');
const WebSocket = require('ws');
const path = require('path');
	var perc=0.,per=0.; //计算pi第几位的进度%
	var fperc=0.,fper;  //计算求和公式的进度% , 用到的时间计数,大于间隔时间则输出
	var epsf=1000;   //时间间隔1秒
	var testbig=0000;
	var u=1;
	var RADIX=3;
	var RAD=(1<<RADIX);
	var PS=10;           //模拟16进制小数点后7位
	var NHX=1;          //可靠精度只有8位,!!!万一出现FFFFFFFF的情形,计算也是不可靠的
	var MAX=(u<<((PS-6)<<2))*(u<<((NHX+6)<<2));  //乘子=16^15
const app = express();
let count = 0; 

let totalSum = 0;
let clients = []; // 用于保存客户端连接


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
    
    // 检查是否达到了80个字符的倍数（即模80等于0）  
    if (count % 80 === 0) {  
        // 如果是，则输出当前累积的result，并重置result和count  
        console.log(result0);  
        result0 = ''; // 重置result，为下一行做准备  
        count = 0; // 重置计数器  
    }  
});

console.log('Server running on http://localhost:3000');
