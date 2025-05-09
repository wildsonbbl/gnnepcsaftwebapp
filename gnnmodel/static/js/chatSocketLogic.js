// Set up WebSocket event handlers
function setupChatSocketHandlers() {
  chatSocket.onopen = function (e) {
    console.log("Connected to chat server");
    // Request list of sessions
    setTimeout(function () {
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
    }, 500); // Small delay to ensure connection is fully established
  };

  chatSocket.onmessage = function (e) {
    var data = JSON.parse(e.data);

    // Handle different types of messages
    if (data.action) {
      handleActionMessage(data);
    } else if (data.text) {
      handleChatMessage(data.text);
    }
  };

  chatSocket.onclose = function (e) {
    console.log("Socket closed unexpectedly, please reload the page.");
  };
}

// Handle action messages (session management)
function handleActionMessage(data) {
  switch (data.action) {
    case "sessions_list":
      populateSessionsList(data.sessions);
      break;
    case "session_loaded":
      // Set the current session when initially loaded
      currentSessionId = data.session_id;
      currentSessionName = data.name;
      document.getElementById("current-session-name").textContent =
        currentSessionName;

      // Handle model information
      if (data.model_name) {
        currentModelName = data.model_name;
        document.getElementById("current-model-name").textContent =
          currentModelName;
      }
      // Populate available models if provided
      if (data.available_models && Array.isArray(data.available_models)) {
        availableModels = data.available_models;
        populateModelsList(availableModels, currentModelName);
      }
      // Update MCP Config Path and Button Title
      if (data.mcp_config_path) {
        mcpConfigPath = data.mcp_config_path;
        const mcpConfigFilePathDisplay = document.getElementById(
          "mcp-config-file-path"
        );
        if (mcpConfigFilePathDisplay) {
          mcpConfigFilePathDisplay.textContent = mcpConfigPath;
        }
      }

      // Populate available tools if provided
      if (data.available_tools && Array.isArray(data.available_tools)) {
        availableTools = [...data.available_tools]; // Use spread syntax for a new array
      }

      // Update selected tools based on the potentially filtered list from the server
      if (data.selected_tools && Array.isArray(data.selected_tools)) {
        selectedTools = [...data.selected_tools];
      } else {
        // Fallback if selected_tools isn't provided, select all available
        selectedTools = [...availableTools];
      }
      populateToolsList(availableTools); // Populate/Repopulate tools list

      if (data.tool_descriptions) {
        Object.assign(toolDescriptions, data.tool_descriptions);
      }
      break;
    case "available_models_updated": // New case
      if (data.available_models && Array.isArray(data.available_models)) {
        availableModels = data.available_models;
        populateModelsList(availableModels, currentModelName);
        if (
          !availableModels.includes(currentModelName) &&
          availableModels.length > 0
        ) {
          // currentModelName is no longer valid, UI will show "Model"
          // User needs to select a new one.
        }
        showToast("Models list updated.", "success");
      }
      break;
    case "ollama_offline": // New case
      showToast(
        "Ollama is offline. Please start Ollama on http://localhost:11434/ and try again.",
        "warning"
      );
      break;
    case "tools_changed":
      if (data.selected_tools) {
        selectedTools = [...data.selected_tools];
      }
      // Update available tools if server sends an updated list
      if (data.available_tools && Array.isArray(data.available_tools)) {
        availableTools = [...data.available_tools];
      }
      populateToolsList(availableTools); // Repopulate with potentially new available tools
      showToast("Tools selection changed successfully");
      break;
    case "mcp_activated":
      availableTools = [...data.available_tools];
      selectedTools = [...data.selected_tools]; // Update selected tools based on server response
      populateToolsList(availableTools);
      if (data.activated_tools && data.activated_tools.length > 0) {
        showToast(
          `MCP Servers activated. New tools available: ${data.activated_tools.join(
            ", "
          )}`
        );
        showToast("To terminate MCP server, click the AI Chat button");
      } else {
        showToast(
          "MCP Servers processed. No new tools were activated (check config or server logs)."
        );
      }
      if (data.tool_descriptions) {
        Object.assign(toolDescriptions, data.tool_descriptions);
      }
      break;
    case "mcp_activation_failed":
      showToast(`MCP Server activation failed: ${data.error}`, "error");
      break;
    case "model_changed":
      // Update the current model when changed
      if (data.model_name) {
        currentModelName = data.model_name;
        document.getElementById("current-model-name").textContent =
          currentModelName;
        // Update the active class in the models dropdown
        updateActiveModel(currentModelName);
      }
      break;
    case "session_created":
      // Close the modal
      var modal = bootstrap.Modal.getInstance(
        document.getElementById("newSessionModal")
      );
      if (modal) modal.hide();
      // Atualiza a lista de sessões
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
      break;
    case "load_messages":
      messages = data.messages;
      updateChatLog();
      break;
    case "session_renamed":
      if (data.session_id === currentSessionId) {
        currentSessionName = data.name;
        document.getElementById("current-session-name").textContent =
          currentSessionName;
      }
      // Refresh sessions list
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
      // Close the modal
      var modal = bootstrap.Modal.getInstance(
        document.getElementById("renameSessionModal")
      );
      if (modal) modal.hide();
      break;
    case "end_turn":
      const generatingContainer = document.getElementById("bottom-chat-log");
      generatingContainer.innerHTML = "";
      generatingContainer.scrollIntoView();
      if (data.type === "stop_action") {
        showToast("Stopped generating", "error");
      }
      if (data.type === "interrupted") {
        showToast("LLM interrupted generating", "error");
      }
      break;
    case "ongoing_turn":
      showGeneratingMessage();
      break;
    case "session_deleted":
      if (data.success) {
        // If the current session was deleted, update UI
        if (data.session_id === currentSessionId) {
          messages = [];
          updateChatLog();
        }

        // Close the modal
        var modal = bootstrap.Modal.getInstance(
          document.getElementById("deleteSessionModal")
        );
        if (modal) modal.hide();

        // Show success message
        showToast("Session deleted successfully");
      } else {
        showToast("Failed to delete session", "error");
      }
      break;
    case "mcp_config_content":
      const mcpConfigContentInput = document.getElementById(
        "mcp-config-content-input"
      );
      const mcpConfigError = document.getElementById("mcp-config-error");
      mcpConfigError.classList.add("d-none");
      mcpConfigError.textContent = "";
      if (data.error) {
        mcpConfigContentInput.value = "";
        mcpConfigError.textContent = `Error loading MCP config: ${data.error}`;
        mcpConfigError.classList.remove("d-none");
      } else {
        mcpConfigContentInput.value = data.content;
      }
      break;
    case "mcp_config_saved":
      if (data.success) {
        showToast("MCP configuration saved successfully.");
        var modal = bootstrap.Modal.getInstance(
          document.getElementById("mcpConfigModal")
        );
        if (modal) modal.hide();
        // Optionally, you might want to re-fetch session data if config changes affect available tools immediately
        // For now, activating MCP will re-fetch.
      } else {
        const mcpConfigError = document.getElementById("mcp-config-error");
        mcpConfigError.textContent = `Error saving MCP config: ${data.error}`;
        mcpConfigError.classList.remove("d-none");
        showToast(`Failed to save MCP configuration: ${data.error}`, "error");
      }
      break;
    case "api_key_required":
      showToast(
        `Google API Key is required for the current model (${data.model_name}). Please set the GOOGLE_API_KEY or GEMINI_API_KEY environment variable.`,
        "error"
      );
      break;
  }
}

// Handle chat messages
function handleChatMessage(message) {
  if ((message.source == "assistant") | (message.source == "user")) {
    messages.push(message);
  }
  updateChatLog();
}

// Add a function to change the model
function changeModel(modelName) {
  if (modelName === currentModelName) return;

  chatSocket.send(
    JSON.stringify({
      action: "change_model",
      model_name: modelName,
      tools: selectedTools,
    })
  );
}
