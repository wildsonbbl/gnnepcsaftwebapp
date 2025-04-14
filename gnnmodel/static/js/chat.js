var wss_protocol = window.location.protocol == "https:" ? "wss://" : "ws://";
var chatSocket = new WebSocket(
  wss_protocol + window.location.host + "/ws/chat/"
);
var messages = [];

chatSocket.onopen = function (e) {
  console.log("Connected to chat server");
};

chatSocket.onmessage = function (e) {
  var data = JSON.parse(e.data);
  var message = data["text"];
  if ((message.source == "assistant") | (message.source == "user")) {
    messages.push(message);
  }

  var str = "";
  messages.forEach(function (msg) {
    str += `<div class="d-flex flex-row mb-4 ${
      msg.source == "assistant"
        ? "justify-content-start"
        : "justify-content-end"
    }">
        <div class="p-3 ms-3  ${
          msg.source == "assistant" ? "bot-message" : "user-message"
        }">
        <p class="small mb-0">${msg.msg}</p></div></div>`;
  });

  str += ` <div class="d-flex flex-row mb-4 justify-content-center" id="bottom-chat-log">
      <p class="small mb-0">${
        message.end_turn == false ? "Generating response..." : ""
      }</p>
      </div>`;

  document.querySelector("#chat-log").innerHTML = str;
  document.getElementById("bottom-chat-log").scrollIntoView();
};

chatSocket.onclose = function (e) {
  console.log("Socket closed unexpectedly, please reload the page.");
};

document.querySelector("#chat-message-input").focus();
document.querySelector("#chat-message-input").onkeyup = function (e) {
  if (e.keyCode === 13) {
    // enter, return
    document.querySelector("#chat-message-submit").click();
  }
};

document.querySelector("#chat-message-submit").onclick = function (e) {
  var messageInputDom = document.querySelector("#chat-message-input");
  var message = messageInputDom.value;
  chatSocket.send(
    JSON.stringify({
      text: message,
    })
  );

  messageInputDom.value = "";
};

document
  .getElementById("chat-log-save_button")
  .addEventListener("click", function () {
    const dataStr =
      "data:text/json;charset=utf-8," +
      encodeURIComponent(JSON.stringify(messages));
    const downloadAnchorNode = document.createElement("a");
    downloadAnchorNode.setAttribute("href", dataStr);
    downloadAnchorNode.setAttribute("download", "messages.json");
    document.body.appendChild(downloadAnchorNode);
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
  });
