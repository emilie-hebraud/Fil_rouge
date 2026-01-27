
# import libraries

def train_step(model: torch.nn.Module,
               dataloader: torch.utils.data.DataLoader,
               loss_fn,
               metric,
               optimizer,
               device):
    
    model.train()
    metric.reset()
    full_score, full_loss, counts = 0,0,0
    
    for x_batch, y_batch in train_loader:
        x_batch = x_batch.to(device)
        y_batch = y_batch.to(device)
        logits = model(x_batch)

        optimizer.zero_grad()
        loss = loss_fn(logits, y_batch.float())
        loss.backward()
        optimizer.step()
        probs = torch.sigmoid(logits)
        preds = (probs > 0.5).int()
        
        # Metrics
        full_loss += loss.item() * x_batch.size(0)
        metric.update(preds.squeeze(1), y_test.int().squeeze(1))
        counts += x_batch.size(0)
    score = metric.compute()
    metric.reset()
    return full_loss, score, counts
    

def test_step(model: torch.nn.Module,
              dataloader: torch.utils.data.DataLoader,
              loss_fn,
              metric,
              device):
    model.eval()
    metric.reset()
    with torch.no_grad():
                full_loss, full_score, counts = 0, 0, 0
                for x_test,y_test in test_loader:
                    
                    print(x_test.shape)
                    print(x_test.size(0))
                    
                    x_test,y_test = x_test.to(device), y_test.to(device)
                    y_logits = model(x_test)
                    probs = torch.sigmoid(y_logits)
                    preds = (probs > 0.5).int()
                    loss = loss_fn(y_logits, y_test.float())
                    full_loss += loss.item() * x_test.size(0)
                    metric.update(preds.squeeze(1), y_test.int().squeeze(1))
                    counts += x_test.size(0)
    score = metric.compute()
    metric.reset()
    return full_loss, score, counts